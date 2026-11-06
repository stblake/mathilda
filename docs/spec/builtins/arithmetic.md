# Arithmetic

Core arithmetic: the basic operators (`Plus`, `Times`, `Power`, `Sqrt`), exact and approximate numbers (`N`, precision control), combinatorial functions (`Factorial`, `Binomial`, `Fibonacci`), digit and radix manipulation, and complex-number component functions. Number-theoretic functions (`GCD`, `PowerMod`, `PrimeQ`, `FactorInteger`, ...) live in [`number-theory.md`](number-theory.md).

## Numeric literals (parser-level)


Numeric literals are classified by the parser based on the precision
the user typed.

- **Integer / BigInt** — literals with no `.`, `e`, `E`, or `*^`.
  Digit-count drives the choice between `EXPR_INTEGER` (`int64_t`) and
  `EXPR_BIGINT` (GMP).
- **`EXPR_REAL`** — real literals of at most **17 significant digits**
  (enough to name every double), counted in the literal text with leading
  zeros excluded, as in the reference: `Precision[1.2887955878011237]` and
  `Precision[9.9999999999999999]` are `MachinePrecision` (v0.6136; the cut
  used to be an implied precision of 15.95 digits, so 17-digit literals
  were arbitrary precision).
- **`EXPR_MPFR`** (arbitrary precision) — real literals of 18 or more
  significant digits (`1.00000000000000001`). The MPFR bit count is
  `ceil(implied_precision × log2(10))` and the literal text is fed
  directly to `mpfr_set_str`, so no precision is lost via a double
  round-trip.

Implied precision is measured from the **mantissa only**, ignoring any
`e` / `E` / `*^` exponent suffix:

    implied_precision = digits_after_decimal + log10(|mantissa_value|)

Consequences:

- `1.0e22` keeps `MachinePrecision` — the mantissa `1.0` has implied
  precision 1, regardless of the large exponent.
- `1.123456789012345678e22` is `EXPR_MPFR` — the mantissa has implied
  precision ≈ 18.
- Pure integer-shaped reals (`1e10`, no `.`) always stay machine
  precision.
- Zero literals (`0.0`, `0.000`, `0.0e10`) always stay machine
  precision regardless of the typed digits.

Examples (verified against the reference):

```mathematica
In[1]:= Precision[-29037945.290347]
Out[1]= MachinePrecision      (* implied ≈ 13.46, below threshold *)

In[2]:= Precision[
            -29037852093587905730945.29034875093457832094573984537498]
Out[2]= 54.4864               (* MPFR @ 181 bits *)

In[3]:= Accuracy[
            -29037852093587905730945.29034875093457832094573984537498]
Out[3]= 32.0235
```

Reference-style backtick suffixes (`` ` ``, `` `` ``) still override
the implicit decision: `3.14`50` is 50-digit MPFR, `3.14``49` is
49-digit accuracy, and `3``50` is `EXPR_REAL`.

## N

Numerical evaluation.
- `N[expr]` — machine-precision approximation.
- `N[expr, p]` — `p`-digit (arbitrary-precision, MPFR) approximation.
- `N[expr, {p, a}]` — precision goal `p` and accuracy goal `a`. The precision
  goal decides, and an `Infinity` precision goal leaves the accuracy goal as the
  digit count. `N` is `Protected` only, not `Listable` (as in the reference), so
  the pair is not threaded over: `N[Pi, {MachinePrecision, MachinePrecision}]`
  is machine `Pi`.

**N-values**: a rule for `N` attached to a symbol or a head is used wherever
that symbol or call occurs, not only at the top: after `z2 /: N[z2] = 3/2` and
`h /: N[h[x_]] := 2 x`, `N[2 + z2]` is `3.5` and `N[1 + h[3]]` is `7.`. The
assignment `N[x] = v` (or `N[f[...]] := v`) stores the same rule on `x` (or on
`f`), since `N` itself is `Protected`. The rules are the symbol's `NValues`
(never its `UpValues`). A rule `N[x] = v` is listed, and matched, as
`N[x, {MachinePrecision, MachinePrecision}] :> v`, as in the reference;
`NValues[x] = {...}` replaces the list and `Clear` removes it. A rule's value is numericalized itself
(`N[x] = 5` makes `N[x]` `5.`). A rule written without a precision argument
serves machine-precision `N` only, as the reference's N-value for
`{MachinePrecision, MachinePrecision}` does.

**Held arguments**: in a position its head holds (`If`'s branches, a
`Function` body, any `HoldFirst`/`HoldRest`/`HoldAll` slot) `N` numericalizes the
structure but does not evaluate it: `N[If[y == 2, Print[a], y + 1]]` is
`If[y == 2., Print[a], y + 1.]` and prints nothing. The same holds for the
numericalization `Equal` does on a symbolic side, so `If[...] == 2` never runs
a held branch.

**Adaptive precision**: `N[expr, p]` aims for `p` good digits. Cancellation
in a sum loses leading digits, so when an addition during the evaluation lost
more bits than the working precision's guard (the part not already planned for
a large exact leaf or exponent), `N` evaluates again at a higher working
precision: at least 64 bits more, and directly to `p` plus the measured loss
when the pass lost fewer bits than it had (after a cancellation to an exact
zero, which measures nothing, the extra precision doubles). A result is
*certified* when (a) it agrees to `p` digits with an evaluation at least 32
bits lower -- directly, or by extrapolation, since a fixed loss of `L` bits
leaves `W - L` good bits at every working precision `W` -- and (b) the worst
cancellation of its own pass left `p` bits. (b) is what keeps an exact `0.`
from certifying another exact `0.`: `N[E^1000 (E^(1/1000 - E^-1000) -
E^(1/1000)), 30]` is a zero at every precision below ~1450 bits and was
answered `0.0`; (a)'s gap is what keeps a short last step to the cap from
agreeing on rounding noise.

The limit is `p + $MaxExtraPrecision` digits (default 50, settable:
`Block[{$MaxExtraPrecision = 1000}, N[expr, 30]]`). At the limit a part that
is certified is kept; a part that is zero, or keeps shrinking as the precision
rises, is `0.` -- a value of size `2^(T - W)` at every `W` cannot be told
from an exact zero times `2^T` by evaluation -- and any other unsettled part is
kept at the digits it can claim; in both cases (and for an `Overflow[]` or
`Indeterminate` at the limit, a division by such a zero) `N` emits
`N::meprec` ("Internal precision limit $MaxExtraPrecision = 50 reached while
evaluating ..."), which `Check` and `Quiet` see. The message for a vanished
part is `N`'s own: the numericalization of an exact operand inside an inexact
`Plus`/`Times` (`N[1, 30] + (Sin[Pi/7]^2 + Cos[Pi/7]^2 - 1)`) keeps it as a
silent `0.`, since only the partner's precision is needed there. If a higher precision comes
back worse (`Overflow[]`, a non-finite value) than a lower one, `N` stops and
keeps the last good result. An expression that cancels nothing is evaluated
once. A power `b^n` with a huge exact `n` (left unevaluated because it is
huge) multiplies the base's rounding error by `n`, so `N[b^n, p]` plans
`log2|n|` extra bits for it. Examples: `N[Cos[10^-20] - 1, 20]` is
`-4.99999999999999999999*10^-41`; `N[(Exp[10^-30] - 1) 10^30 - 1, 20]` is
`5.*10^-31` (two 100-bit cancellations); `N[E^1000 (E^(1/1000 - E^-1000) -
E^(1/1000)), 30]` needs ~435 extra digits, so it is `0.` with `N::meprec`, and
`-1.00100050016670834166805575399` under `$MaxExtraPrecision = 1000`;
`N[(Log[Cos[x]] + x^2/2 + x^4/12 + x^6/45)/x^8 /. x -> 10^-10, 30]` needs ~82
extra digits (was 18 good digits reported as 30, silently);
`N[Sin[Pi/7]^2 + Cos[Pi/7]^2 - 1, 40]` is `0.` with `N::meprec`.
Machine-number inputs gain no digits, so `N[0.1 + 0.2 - 0.3, 32]` stays
`5.55112*10^-17`. Only the `N` builtin's `N[expr, p]` adapts; internal
working-precision numericalization evaluates once. Arithmetic on numbers that
are already inexact does not track cancellation: with `a = N[E^1000, 30]`,
`a - (a - 1)` is `0.` and its `Precision` is the operands' 30.1.

**Precision semantics**:
- `Precision[expr]` and `Accuracy[expr]` are one number for any expression,
  the **minimum** over its numbers (they are not `Listable`): exact parts
  impose nothing, a machine number makes it `MachinePrecision`.
  `Precision[N[HilbertMatrix[4], 30]]` is `30.`, `Precision[{1, 2.5`20, N[1, 30]}]`
  is `20.`, `Precision[{1, 2.5, N[1, 30]}]` is `MachinePrecision`,
  `Precision[{1, 2}]` and `Precision[x]` are `Infinity`,
  `Accuracy[{1.5, 100.`30}]` is `15.7785`. A packed array is read from its
  buffer (an integer array is exact: `Infinity`); a `SparseArray` takes the
  minimum over its background and stored values. (Through v0.6225 both were
  `Listable`, so a matrix answered a matrix of precisions.)
- **Machine rounding is IEEE round-to-nearest, ties to even, with gradual
  underflow** (v0.6179). An exact rational becomes the double a correctly
  rounded decimal literal of it would give, through exactly one rounding
  (`numeric_mpq_to_double`, in GMP integers). The earlier path rounded twice
  (to a 53-bit MPFR, then to the double), which could be off by one ulp below
  `$MinMachineNumber`, where a subnormal holds fewer than 53 bits:
  `N[3/2 2^-1074 - 2^-1200]` gave `2^-1073` instead of `2^-1074`. An `int64`
  numerator or denominator past `2^53` also no longer rounds before the
  division: `N[(2^62 + 1)/3]` is the nearest double. `Rationalize[x, 0]` and
  `SetPrecision[x, Infinity]` therefore round-trip through `N` exactly.
- `SetPrecision[x, Infinity]` (and `SetAccuracy[x, Infinity]`) replaces each
  approximate number by its **exact** value. A machine or MPFR number is a
  dyadic rational, so nothing rounds: `SetPrecision[0.1, Infinity]` is
  `3602879701896397/36028797018963968`, `SetPrecision[5.*^-324, Infinity]` is
  `2^-1074`. Exact numbers, symbolic constants and symbols are unchanged
  (`SetPrecision[{0.5, Pi, x + 0.25}, Infinity]` is `{1/2, Pi, 1/4 + x}`), and
  held forms are not entered. (Until v0.6179 it acted as `MachinePrecision`,
  which numericalised `Pi` and left `0.1` inexact; `SetAccuracy[x, Infinity]`
  did not evaluate.)
- `N[expr]` gives a **machine-precision** approximation of the whole result. It
  targets machine precision even when `expr` is already an approximate
  arbitrary-precision number, so `N[N[Pi, 100]]` (i.e. `N[Pi, 100] // N`)
  collapses back to `MachinePrecision` rather than keeping 100 digits, matching
  the reference. Likewise `N[2.5`100]` comes back at `MachinePrecision`. To
  request a specific precision, use the two-argument form.
- `N[expr, p]` is an explicit precision request, but like `N[expr]` it **never
  increases** the precision of a number that is already approximate — `N`
  cannot manufacture digits that aren't there. Exact parts (integers,
  rationals, constants) are produced at `p` digits, while each inexact leaf is
  capped at `min(existing, p)`: a machine real stays `MachinePrecision`
  (`N[1., 50]` → machine), a `1.0`30` stays 30 digits under `N[.., 50]`, and an
  already-higher-precision value is still **lowered** (`N[N[Pi, 100], 20]` → 20
  digits). Contrast `SetPrecision[x, p]`, which *does* pad an approximate value
  up to `p` digits. Example: `N[{1., 1, 1.0`30}, 50]` has precisions
  `{MachinePrecision, 50.272, 30.103}`.
- Inexact contagion is unaffected: mixing a machine real with a
  higher-precision value collapses to machine precision
  (`1. + N[Pi, 100]` → `4.14159`), since `MachinePrecision` is the floor.
- Contagion is **narrower than explicit `N`**: an approximate operand in a
  `Plus`/`Times` numericalizes only those siblings that are *numeric*
  quantities (`1. + Sin[2]` → `1.9093`, `1. + Pi + x` → `4.14159 + x`), as
  in the reference (v0.6019). A sibling that is not numeric keeps its exact
  parts: `1. + 2 x` stays `1. + 2 x` (it was `1. + 2. x`), `1. - x` is not
  `1. - 1. x`, `1. + x - x` is `1.` (it was `1. + 0. x`), and `1.5 (2 + a)`
  stays as it is. Nor does contagion thread into the
  arguments of a non-numeric head. So `0.333 x[1]` stays `0.333 x[1]` (the `1`
  is `x`'s argument, not an arithmetic operand) and `1. Abs[x[3]]` keeps
  `x[3]`. Explicit `N[]` still threads everywhere, so `N[x[1]]` is `x[1.]`.
- Rationals whose numerator or denominator overflow `int64`/`double` still
  numericalize to a **single** real, not a frozen `Rational[Real, Real]`
  (`N[1/10^30]` → `1.e-30`, `N[10^400/3]` → `3.33e+399`). Out-of-range machine
  results are carried as a machine-precision MPFR real so the arbitrary
  exponent survives (matching `N[1001!]`).

**Working precision** — `N` does not implement numeric math. It replaces
*leaves* with numeric values and hands the result back to the evaluator, so the
leaf conversion is the one step that must never be what loses information.
Rounding an exact leaf to a double costs a relative 2⁻⁵³, and the function then
amplifies that by its condition number — |x| for the trig family, which for
`Sin[3141592653589793238]` is enough to destroy every digit. `N` therefore picks
a working precision from the input, evaluates there, and rounds **once** at the
end. Both terms are derived, not tuned:

- an exact integer (or dyadic rational) is carried at its own mantissa width, at
  which it is represented exactly and contributes no error at all;
- every other exact rational is carried at `p + ⌈log₂|x|⌉ + 64` bits, since a
  relative 2⁻ᵂ on an argument of size |x| reaches the answer amplified by |x|.

Values that already fit a double exactly — `10^20`, `2^53`, `3/4`, `1/3` — keep
taking the plain machine path unchanged. The guarantee is that `N[f[x]]` for
exact `x` is accurate to the precision it reports, for every builtin:

```mathematica
In[1]:= N[Sin[3141592653589793238]]
Out[1]= -0.446315

In[2]:= N[Sin[3141592653589793238], 30]
Out[2]= -0.4463151633593201122016036193238
```

An *inexact* argument is a different question and is answered differently: a
machine real is its own exact binary value, so `Sin[1.0*^25]` is the sine of
that double (`-0.305258`), not of 10²⁵ (`-0.744790`).

**Machine numbers have an arbitrary exponent.** A machine number carries a
53-bit mantissa but is not restricted to the IEEE exponent range, so
`N[Exp[1000]]` is `1.97007×10⁴³⁴`, not `Infinity`, and `N[Exp[-1000]]` is
`5.07596×10⁻⁴³⁵`, not `0.` — the same representation `N[1001!]` has always
used. Past MPFR's own exponent range (~10³²³²²⁸⁴⁵⁸) the result is genuinely
infinite and is reported as such.

```mathematica
In[1]:= N[Pi, 100] // N
Out[1]= 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068

In[2]:= Precision[%]
Out[2]= 100.243
```

## Plus (+)

Symbolic sum.
- `a + b + ...`

**Features**:
- `Flat`, `Orderless`, `Listable`.
- Combines numeric constants and collects like terms (e.g., `x + 2x` becomes `3x`).
- Like-term collection is hash-based (O(N) expected) and canonically sorts only
  the collapsed result, so a large flat sum (e.g. `Total` of 10^5 monomials) does
  not pay the generic `Orderless` sort over all N raw terms. `Times` collects
  factors the same way.
- Honours a cleared `Orderless` (`ClearAttributes[Plus, Orderless]`), as the
  reference does: numbers still collect into one leading term, but like terms
  combine only with the term just before them and the result keeps argument
  order -- `b + a` stays `b + a`, `b + a + 1 + 2` is `3 + b + a`,
  `b + a + b` stays, `b + a + a` is `b + 2 a`. `Times` likewise (`b a b`
  stays, `b a^2 a` is `b a^3`, `b a 2 3` is `6 b a`). With the default
  attributes the behaviour and cost are unchanged.
- A collected term is a flat `Times`: `2 a b + x` holds `Times[2, a, b]`, never
  `Times[2, Times[a, b]]`. The coefficient is spliced in front of the
  monomial's factors.
- Returns `0` if no arguments are provided.
- Returns `Overflow[]` if integer addition overflows or if any argument is `Overflow[]`.
- A finite term plus a real `±Infinity` absorbs the finite part (`Infinity + 3 → Infinity`),
  but a finite term plus a **directed (complex) infinity** keeps its direction rather
  than collapsing to a real `Infinity`: `1/2 - I Infinity` stays
  `Plus[1/2, Times[Complex[0, -1], Infinity]]` (so a Mellin–Barnes / Bromwich
  vertical line `s = c ± I Infinity` remains recognisable to `Integrate`).

```mathematica
In[1]:= 1 + 2 + x + 2*x
Out[1]= 3 + 3 x

In[2]:= a + b - (a + b)
Out[2]= 0

In[3]:= Interval[{1, 2}] + Interval[{-1, 3}]
Out[3]= Interval[{0, 5}]

In[4]:= NDArray[{1., 2., 3.}] + 10
Out[4]= NDArray[{11.0, 12.0, 13.0}]
```

## Total

Gives the total of elements in a list.
- `Total[list]`
- `Total[list, n]`
- `Total[list, {n}]`
- `Total[list, {n1, n2}]`

**Features**:
- `Protected`.
- `Total[list]` is equivalent to `Apply[Plus, list]`. The list elements are
  already evaluated, so the level-1 sum collects like terms directly (skipping a
  redundant re-evaluation and sort of every term); summing a large symbolic list
  is fast even when the result collapses to a few distinct terms.
- `Total[list, n]` totals all elements down to level `n`.
- `Total[list, {n}]` totals elements at level `n` only.
- Supports negative levels to count from the bottom (`-1` is the last dimension).
- Handles ragged arrays correctly by summing from the inside out when multiple levels are specified.
- `Total[list, Infinity]` totals all atoms in the expression.
- `Total[{}]` is `0` — the additive identity, at any level spec, matching
  `Plus @@ {}`. It answered `{}` until 2026-07-31, which then propagated as a
  non-number through anything consuming it (a `Select` that matched nothing is
  the usual way to reach it).
- Reads a packed buffer directly (`ndred_total`), exactly on an `int64` buffer:
  the accumulation abandons the whole result on overflow so the List path
  re-runs it and GMP answers.

```mathematica
In[1]:= Total[{a, b, c, d}]
Out[1]= a + b + c + d

In[2]:= Total[{{1, 2}, {3, 4}}]
Out[2]= {4, 6}

In[3]:= Total[{{1, 2}, {3, 4}}, 2]
Out[3]= 10

In[4]:= Total[{{1, 2}, {3}}, 2]
Out[4]= 6
```

## Accumulate

Gives the running cumulative totals of the elements in a list.
- `Accumulate[list]`
- `Accumulate[list, Method -> "CompensatedSummation"]`

**Features**:
- `Protected`.
- `Accumulate[list]` has the same length as `list`, and is effectively equivalent to `FoldList[Plus, list]`.
- The head of the input is preserved, so `Accumulate[f[a, b, c]]` returns `f[a, a + b, a + b + c]`.
- Threads naturally over rows via `Listable` `Plus`, so `Accumulate` of a matrix accumulates within columns.
- Works on machine integers, GMP arbitrary-precision integers, machine-precision doubles, and symbolic expressions.
- With `Method -> "CompensatedSummation"`, Kahan compensated summation is used in double precision when every element is a machine number, reducing floating-point round-off error. For symbolic or mixed input the option is silently ignored and the standard symbolic accumulation is returned.

```mathematica
In[1]:= Accumulate[{a, b, c, d}]
Out[1]= {a, a + b, a + b + c, a + b + c + d}

In[2]:= Accumulate[{{a, b}, {c, d}, {e, f}}]
Out[2]= {{a, b}, {a + c, b + d}, {a + c + e, b + d + f}}

In[3]:= Accumulate[f[a, b, c, d]]
Out[3]= f[a, a + b, a + b + c, a + b + c + d]

In[4]:= Accumulate[{1, 2, 3, 4, 5}]
Out[4]= {1, 3, 6, 10, 15}

In[5]:= Accumulate[{1.0, 2.0, 3.0}, Method -> "CompensatedSummation"]
Out[5]= {1., 3., 6.}
```

## Differences

Gives the successive differences of the elements in a list. The inverse of `Accumulate`/`FoldList[Plus, ...]`.
- `Differences[list]`
- `Differences[list, n]`
- `Differences[list, n, s]`
- `Differences[list, {n1, n2, ...}]`

**Features**:
- `Protected`.
- `Differences[list]` gives `{list[[2]] - list[[1]], list[[3]] - list[[2]], ...}`, of length `l - 1`.
- `Differences[list, n]` applies the first-difference operator `n` times, giving length `l - n`. `n` must be a non-negative integer (`n = 0` returns `list` unchanged).
- `Differences[list, n, s]` takes differences of elements step `s` apart, of length `l - n |s|`. The step `s` is a nonzero integer; for `s < 0` the elements are subtracted in the opposite order.
- `Differences[list, {n1, n2, ...}]` gives the successive `nk`-th differences at level `k` of a nested list, and is equivalent to `Differences[Differences[list, n1], {0, n2, ...}]`. Each `nk` must be a non-negative integer.
- Subtraction threads element-wise over sublists via the `Listable` `Plus`/`Times`, so for a matrix `m`, `Differences[m]` (= `Differences[m, 1]` = `Differences[m, {1, 0}]`) differences successive rows within each column, while `Differences[m, {0, 1}]` differences columns within each row.
- The head of the input is preserved. A list shorter than the reduction yields the empty list.
- Works on machine integers, GMP arbitrary-precision integers, machine-precision doubles, and symbolic expressions.

```mathematica
In[1]:= Differences[{a, b, c, d, e}]
Out[1]= {-a + b, -b + c, -c + d, -d + e}

In[2]:= Differences[{a, b, c, d, e}, 2]
Out[2]= {a - 2 b + c, b - 2 c + d, c - 2 d + e}

In[3]:= Differences[{1, 2, 4, 8, 16, 32}, 1, 2]
Out[3]= {3, 6, 12, 24}

In[4]:= Differences[{{a11, a12, a13}, {a21, a22, a23}}, {0, 1}]
Out[4]= {{-a11 + a12, -a12 + a13}, {-a21 + a22, -a22 + a23}}

In[5]:= FoldList[Plus, a, Differences[{a, b, c, d, e}]]
Out[5]= {a, b, c, d, e}
```

## Ratios

Gives the successive ratios of the elements in a list — the multiplicative analog of `Differences`. The inverse of `FoldList[Times, ...]`.
- `Ratios[list]`
- `Ratios[list, n]`
- `Ratios[list, {n1, n2, ...}]`

**Features**:
- `Protected`.
- `Ratios[list]` divides successive elements by preceding ones, giving `{list[[2]]/list[[1]], list[[3]]/list[[2]], ...}` of length `l - 1`; `Ratios[{x1, x2}]` gives `{x2/x1}`.
- `Ratios[list, n]` applies the ratio operator `n` times, giving length `l - n`. `n` must be a non-negative integer (`n = 0` returns `list` unchanged).
- `Ratios[list, {n1, n2, ...}]` gives the successive `nk`-th ratios at level `k` of a nested list, and is equivalent to `Ratios[Ratios[list, n1], {0, n2, ...}]`. Each `nk` must be a non-negative integer.
- Division threads element-wise over sublists via the `Listable` `Power`/`Times`, so for a matrix `m`, `Ratios[m]` (= `Ratios[m, 1]` = `Ratios[m, {1, 0}]`) takes ratios of successive rows within each column, while `Ratios[m, {0, 1}]` takes ratios of columns within each row.
- The head of the input is preserved. A list shorter than the reduction yields the empty list. First ratios are constant for a geometric sequence.
- Works on machine integers, GMP arbitrary-precision integers (exact `Rational` results), machine-precision doubles, and symbolic expressions.

```mathematica
In[1]:= Ratios[{a, b, c, d, e}]
Out[1]= {b/a, c/b, d/c, e/d}

In[2]:= Ratios[{a, b, c, d, e}, 2]
Out[2]= {(a c)/b^2, (b d)/c^2, (c e)/d^2}

In[3]:= Ratios[Table[2^i, {i, 0, 10}]]
Out[3]= {2, 2, 2, 2, 2, 2, 2, 2, 2, 2}

In[4]:= Ratios[{{a11, a12, a13}, {a21, a22, a23}}, {0, 1}]
Out[4]= {{a12/a11, a13/a12}, {a22/a21, a23/a22}}

In[5]:= FoldList[Times, a, Ratios[{a, b, c, d, e}]]
Out[5]= {a, b, c, d, e}
```

> **Packed arrays.** `Ratios[list]` divides in place on a real buffer. An
> **integer** buffer takes the ordinary path: `Ratios[{1, 2, 3}]` is
> `{2, 3/2}`, exact `Rational`s that no buffer holds.

## ListGradient

The numerical gradient of a sampled array by finite differences — a
from-first-principles port of `numpy.gradient`. Second-order central differences
in the interior, one-sided at the edges, per axis, for an array of any rank.
- `ListGradient[f]`
- `ListGradient[f, spacing]`
- `ListGradient[f, spacing, opts...]`

**Features**:
- `Protected`.
- Return shape follows numpy: a rank-1 `f` gives a single vector; a rank-`k` `f`
  gives `{g1, ..., gk}`, one same-shape array per axis (ascending axis order); a
  single requested axis gives one array, not a length-1 list.
- Default (`Method -> "Centered"`, `DifferenceOrder -> 2`) reproduces
  `numpy.gradient` exactly: an interior point is a second-order central
  difference and each extreme endpoint a first-order one-sided difference
  (numpy `edge_order = 1` — the edge accuracy is one order below the interior).
- **spacing** (numpy's `varargs`, as one positional argument):
  - absent → unit spacing on every axis;
  - a scalar (number or symbol) → uniform spacing on every axis;
  - a list of one entry per computed axis, each a scalar spacing or a coordinate
    vector of that axis's length (non-uniform, exact Fornberg weights). For a
    vector, a coordinate list of the same length as `f` is used directly.
- **Options**:
  - `Method -> "Centered"` (default), `"Forward"`, or `"Backward"` — the stencil
    placement. Forward/Backward stay order `p` at every point (shifting inward
    near the far edge); Centered drops the two extreme endpoints to order `p-1`.
  - `DifferenceOrder -> p` (accuracy order, default `2`).
  - `WindowLength -> m` (stencil size; default `Automatic` = `p + 1`). It is an
    alias for the order: an `m`-point stencil is order `m - 1`.
  - `Axis -> All` (default) `| a | {a1, ...}` — restrict to particular axes
    (1-based; negatives count from the end).
- The single kernel is Fornberg's finite-difference-weight recurrence, exact on
  non-uniform grids, so integer input gives exact `Rational`s and symbolic input
  a symbolic result. A machine-`Real` array uses the packed/NDArray fast path,
  and `ListGradient[v]` compiles (rank-1 array form) and auto-compiles.

```mathematica
In[1]:= ListGradient[{1, 4, 9, 16, 25}]
Out[1]= {3, 4, 6, 8, 9}

In[2]:= ListGradient[{a, b, c, d, e}]
Out[2]= {-a + b, -1/2 a + 1/2 c, -1/2 b + 1/2 d, -1/2 c + 1/2 e, -d + e}

In[3]:= ListGradient[{1, 4, 9, 16, 25}, {0, 1, 3, 6, 10}]
Out[3]= {3, 17/6, 73/30, 193/84, 9/4}

In[4]:= ListGradient[{{1, 2, 6}, {3, 4, 5}}, Axis -> 2]
Out[4]= {{1, 5/2, 4}, {1, 1, 1}}

In[5]:= ListGradient[{0, 1, 8, 27, 64, 125, 216}, DifferenceOrder -> 4]
Out[5]= {0, 3, 12, 27, 48, 75, 108}
```

> **Packed arrays.** `ListGradient` works straight on a `float64`/`float32`
> **or complex** (`complex64`/`complex32`) buffer (packed or visible `NDArray`),
> one axis at a time via a strided axpy-style pass — the finite-difference
> weights are real, so a complex element is just two contiguous doubles carried
> through the same pass. An **integer** buffer takes the ordinary path: the
> central difference of exact integers is a `Rational` (e.g.
> `(f[i+1]-f[i-1])/2`) that no buffer holds, so it materialises and the exact
> `List` path answers. Arbitrary-precision (`MPFR`) and exact
> (`Integer`/`Rational`) arrays likewise compute on the exact/symbolic `List`
> path, retaining full precision.

## Times (*)

Symbolic product.
- `a * b * ...`

**Features**:
- `Flat`, `Orderless`, `Listable`.
- Combines numeric constants and groups identical bases into `Power` expressions.
- **Numeric bases under one symbolic exponent merge** (v0.594), as in the
  reference: `2^a 3^a` is `6^a`, `3^-a 5^a` is `(5/3)^a`, `2^a 3^a 5^-a` is
  `(6/5)^a`, and `Sinh[1]^a Cosh[1]^-a` is `Tanh[1]^a`. Since v0.6016 any
  exponent that is not an integer, a rational or a machine real counts, so a
  numeric one merges too: `2^Sqrt[2] 3^-Sqrt[2]` is `(2/3)^Sqrt[2]`, `2^Pi 3^Pi`
  is `6^Pi`, `2^I 3^-I` is `(2/3)^I` (rational exponents are the radical
  fusion's: `Sqrt[2]/Sqrt[3]` is `Sqrt[2/3]`). The exponents must be
  equal or exactly opposite, and every base a numeric quantity. Symbolic bases
  never merge: `x^a y^a` and `2^a f[1]^a` stay. `(b1 b2)^e = b1^e b2^e` needs
  all but one base to be a positive real, so a group admits at most one other
  base: `(-2)^a 3^a` is `(-6)^a`, while `(-2)^a (-3)^a` stays.
  - The merged power keeps the exponent of its non-positive base, if there is
    one, because only a positive base may be inverted to flip a sign.
  - Otherwise a group whose exponents all agree keeps that exponent
    (`3^-a 5^-a` is `15^(-a)`), and a mixed group takes the sign without a
    negative leading coefficient (`a`, not `-a`). A sum takes the sign whose
    last term is negative (v0.6016): `2^(1 + a) 3^(-1 - a)` is
    `(3/2)^(-1 - a)`, `2^(a - 1) 3^(1 - a)` is `(3/2)^(1 - a)`, and
    `2^(1 + Sqrt[3]) 3^(-1 - Sqrt[3])` is `(3/2)^(-1 - Sqrt[3])`. A last term
    with a complex coefficient reads as negative when the printer writes it
    with a minus (`-I t`, `(-1 + I) t`); a product flips only for a negative
    imaginary coefficient (`2^(-I t) 3^(I t)` is `(3/2)^(I t)`, while
    `(2/3)^((-1 + I) t)` stays); a complex number keeps the sign with a
    positive imaginary part (`2^(-I) 3^I` is `(3/2)^I`). Exponents with an
    inexact number do not merge (`2^(2.5 t) 3^(-2.5 t)` stays).
- **A power of a product splits off its positive numeric factors** (v0.6005),
  as in the reference, when it is a factor of a larger product and its
  exponent is not a number:
  - `y (2 x)^b` is `2^b x^b y`, and `x^a (Pi x)^b` is `Pi^b x^(a + b)`.
  - `(x/2)^b` contributes `2^-b`; `(3/2 x)^b` keeps `(3/2)^b`.
  - A negative number contributes its absolute value: `(-2 x)^b` gives
    `2^b (-x)^b`. A pure-imaginary coefficient does the same: `(2 I x)^b`
    gives `2^b (I x)^b`.
  - The rest of the base stays as one power: `(2 x y)^b` gives `2^b (x y)^b`.
  - A base that is all numeric, such as `(2 Pi)^b`, splits only when another
    factor shares one of its bases. So `Pi^a (2 Pi)^b` is `2^b Pi^(a + b)`
    and `Pi^(D/2) (2 Pi)^-D` is `2^-D Pi^(-D/2)`, while `y (2 Pi)^b` and
    `(-2 Pi)^b y` stay.
  - A lone `(2 x)^b` is not a product and stays.
  - The split factors are spliced into one flat product, so like bases
    combine before the numeric-base merge above can rebuild `(2 Pi)^b`.
- **Numeric powers whose bases share a factor are written over one coprime
  basis** (v0.6008) when an exponent is symbolic, as in the reference.
  - Examples: `4^a 6^b` is `2^(2 a + b) 3^b`, `8 2^a 4^b` is
    `2^(3 + a + 2 b)`, `12 2^n` is `3 2^(2 + n)`, `4^a/2` is `2^(-1 + 2 a)`, and
    `(2/3)^a (9/4)^b` is `(2/3)^(a - 2 b)`.
  - Powers of one base add their exponents first (`2 2^a` is `2^(1 + a)`),
    then powers with equal exponents are multiplied, so `2^a 8^a` is `16^a`;
    the basis is taken of the values that remain (v0.6225), so `6^x 4^x 9` is
    `8^x 3^(2 + x)` (24 and 9 refine to {3, 8}) and `2^a 8^a 4^b` is
    `4^(2 a + b)`. The equal-exponent merge then recombines what ends up equal
    or opposite.
  - **Negative rational bases** (v0.6225, issue #126) join through their
    magnitude: `(-n)^e` is `(-1)^e n^e` for every complex `e`. The sign goes
    back onto a basis element of the power's own value whose exponent is
    still `e`, else into one `(-1)^(e1 + e2 + ...)`: `(-2)^x/2` is
    `(-1)^x 2^(-1 + x)`, `(-6)^x/3` is `(-2)^x 3^(-1 + x)`, `(-6)^x/4` is
    `(-3)^x 2^(-2 + x)`, `(-2)^x (-8)^y` is `(-1)^(x + y) 2^(x + 3 y)`, and
    `(-9)^(1 + x) (-12)^x 4` is `(-1)^(1 + 2 x) 3^(2 + 3 x) 4^(1 + x)`. So
    `(-2)^x/2 - (-1)^x 2^(x - 1)` is `0`.
    - Two negative powers under one exponent never multiply their
      magnitudes: `(-2)^x (-3)^x` stays, `(-2)^x (-4)^x` is
      `(-1)^(2 x) 2^(3 x)`. A negative power's magnitude meets a positive
      power of the same value only to cancel it: `(-1/2)^n 2^n` is `(-1)^n`,
      while `(-2)^x 2^x` is `(-4)^x`.
    - The number factors are multiplied first, and that coefficient merges
      with a negative power of the same base only when that power is the
      first (smallest base) to share a prime with it: `(-3)^x (-3)` is
      `(-3)^(1 + x)` and `(-4)^(2 + n) (-4) 2^(2 x)` is
      `(-1)^(3 + n) 2^(6 + 2 n + 2 x)`, while `(-1/2)^(2 x) (-8)^x (-1/2)` is
      `-(-1)^(3 x) 2^(-1 + x)` and `(-1/2)^x (-1/2) 4` (coefficient `-2`) is
      `-(-1)^x 2^(1 - x)`.
    - Not yet: a negative base under a rational exponent (`(-2)^(1/3)/2`,
      which the reference writes `(-1)^(1/3) 2^(-2/3)`), and Gaussian-integer
      coefficients (`(-2)^x/(2 I)`).
  - The basis comes from gcds alone. A product of huge prime powers such as
    `m^-7 (p/m^4)^(-1 - 3 e) (r/m^3)^(-1 + 4 e)` collapses with no integer
    factoring.
  - Bases that share nothing stay as they are (`3 4^a`, `2^a 3^b`).
  - A shared factor calls for the basis only when at least one of the two
    powers has a symbolic exponent (v0.6016): `17 Sqrt[17]` stays. Rewritten
    as `17^(3/2)` it evaluated straight back, an endless loop through `Times`
    (`$RecursionLimit::reclim`).
- Handles `I` as `Complex[0, 1]`.
- Returns `1` if no arguments are provided.
- Returns `Overflow[]` if integer multiplication overflows or if any argument is `Overflow[]`.
- **Sign distribution.** A product of exactly two factors, a coefficient `-1`
  or `-1/n` and a sum, distributes the sign over the sum, as in the reference:
  `-(-1 + x)` is `1 - x`, `-(x + y) z` is `(-x - y) z` (the parser reads it as
  `Times[Times[-1, x + y], z]`), and `Times[-1/2, x + y]` is `(-x - y)/2`. Any
  other shape stays as it is: `Times[-1, x + y, z]`, `Times[-3/2, x + y]`,
  `Times[-1., x + y]` and `-(x + y)^2`.
- **Sqrt-coefficient absorption**: a non-trivial rational/integer coefficient
  combined with a single `Power[r, -1/2]` group folds into the canonical
  reference form `sign(c) * Sqrt[c^2 / r]`, then Power's existing
  rational-base extraction pulls out any newly-exposed perfect square:
  `14/Sqrt[10] → 7 Sqrt[2/5]`, `78/Sqrt[66] → 13 Sqrt[6/11]`,
  `(3/5)/Sqrt[2/5] → 3/Sqrt[10]`. Skipped when no square is actually
  extracted (e.g. `2/Sqrt[30]` stays as-is to preserve like-term
  collection in `Plus`).

```mathematica
In[1]:= 2 * x * 3 * y
Out[1]= 6*x*y

In[2]:= I * I
Out[2]= -1

In[3]:= 14/Sqrt[10]
Out[3]= 7 Sqrt[2/5]
```

**v0.6035.** A symbolic factor does not rescue `0 * Infinity`: `0 c ComplexInfinity` is `Indeterminate`, with `Infinity::indet`, as `0 ComplexInfinity` is.

## NonCommutativeMultiply (**)

A general associative, but non-commutative, form of multiplication. Use it as a
generalization of ordinary multiplication for special mathematical objects
(e.g. composition in an algebra of differential operators, matrices of
operators, quaternions).
- `a ** b ** c` or `NonCommutativeMultiply[a, b, c]`.

**Features**:
- Attributes `Flat`, `OneIdentity`, `Protected`.
- Instances are automatically flattened (`a ** (b ** c)` and `(a ** b) ** c`
  both become `NonCommutativeMultiply[a, b, c]`), **but no other simplification
  is performed**: arguments are not reordered (`a ** b` stays distinct from
  `b ** a`), `0 ** a` and `1 ** a` are left as given, and `Expand` / `Simplify`
  / `FullSimplify` do not operate on expressions containing it.
- `NonCommutativeMultiply[a]` (one argument) stays unevaluated — unlike
  `Times[a]`, which reduces to `a`.
- Precedence binds tighter than `*` (`Times`) and `.` (`Dot`), and looser than
  `^` (`Power`): `a ** b ^ c` is `a ** (b^c)`, `a b ** c` is `a (b ** c)`.
- No new numeric fast path: `**` is a purely symbolic/structural head, so it has
  no packed/NDArray or `Compile[]` lowering.

```mathematica
In[1]:= a ** b ** c
Out[1]= a ** b ** c

In[2]:= a ** (b ** c) == (a ** b) ** c
Out[2]= True

In[3]:= {a*b == b*a, a**b == b**a}
Out[3]= {True, a ** b == b ** a}

In[4]:= {0 ** a, 1 ** a}
Out[4]= {0 ** a, 1 ** a}

In[5]:= x ** y ** z /. y ** z -> w     (* a rule fires on any contiguous run *)
Out[5]= x ** w

In[6]:= op /: op[i_] ** op[j_] := op[i + j]; x ** op[1] ** op[2]
Out[6]= x ** op[3]
```

Behaviour is added with rules (UpValues on the operands are the usual way). Being
`Flat`, a rule on `**` matches a contiguous, in-order run of operands anywhere in a
longer product, leftmost first; see *Rules on Flat heads* in
[assignment-and-rules.md](assignment-and-rules.md).

## Internal`SyntacticNegativeQ
Does an expression carry a minus sign in the leading position of its printed
form? This is a syntactic test, not a sign test.
- `Internal`SyntacticNegativeQ[expr]`

**Features**:
- True for a negative number or `Rational`, and for a `Complex` whose
  imaginary part is negative. `Complex[1, -2]` and `-I` are `True`;
  `Complex[-1, 2]` is `False`.
- A `Times` or `Plus` is read by its first argument in canonical order:
  `-2 x y`, `-a - b` and `-1 + x` are `True`, while `x - y` and `1 - 2 x` are
  `False` (v0.6012, as in the reference; a sum used to be always `False`).

## Minus
Arithmetic negation as a function.
- `Minus[x]` is `-x`, i.e. `Times[-1, x]`.

**Features**:
- `Listable`, `NumericFunction`, `Protected`, as in the reference.
- Rewrites to `Times[-1, x]`, so numbers, rationals, complex numbers, `Plus`
  distribution and symbolic arguments all behave exactly like `-x`.
- Usable as a sort key: `SortBy[list, Minus]`, `KeySortBy[a, Minus]`.
- Packed arrays and visible `NDArray`s stay on the buffer (`Minus` is on the
  packed-aware list and hands the buffer to `Times`); `Compile[]` lowers it at
  scalar and array shapes.
- Any other argument count emits `Minus::argx` and stays unevaluated.

```mathematica
In[1]:= Minus[3]
Out[1]= -3

In[2]:= Minus[{1, 2.5, 1/2}]
Out[2]= {-1, -2.5, -1/2}

In[3]:= SortBy[{3, 1, 2}, Minus]
Out[3]= {3, 2, 1}
```

## Power (^)

Exponentiation.
- `base ^ exp`

**Features**:
- `Listable`.
- Simplifies integer powers of integers.
- `E^(I q Pi)` with q an integer or half-integer evaluates to ±1 or ±I (`E^(I Pi) = -1`, `E^(-I Pi/2) = -I`), as `Exp[I q Pi]` already did and the reference does. Other rational q stay symbolic (`E^(I Pi/4)`) (v0.458).
- Returns `Overflow[]` if the result exceeds 64-bit integer limits.
- Reduces radicals (e.g., `8^(1/2)` becomes `2*Sqrt[2]`).
- A negative exact coefficient (other than `-1`) in a product under a rational
  power comes out as its magnitude, the sign staying inside:
  `Sqrt[-4 (1 - x^2)]` is `2 Sqrt[-1 + x^2]` and `(-8 a)^(1/3)` is
  `2 (-a)^(1/3)`, as in the reference (`(c w)^p = (-c)^p (-w)^p` for `c < 0`).
  Gated, like the positive case, on the magnitude reducing cleanly
  (v0.591). The reference also splits a coefficient that does not reduce
  (`Sqrt[12 a]` is `2 Sqrt[3] Sqrt[a]`); Mathilda keeps `Sqrt[12 a]`.
- Supports complex results for negative bases (e.g., `(-1)^(1/2)` becomes `I`).
  Higher-power cases for `q == 2` now also reduce: `(-1)^(3/2) → -I`,
  `(-1)^(5/2) → I`, `(-12)^(3/2) → -24 I Sqrt[3]` (the principal-branch
  rule `(-n)^(p/2) = I^p · |n|^(p/2)`).
- For a negative base, the residual `(-1)^(b/q)` exponent is reduced into
  `[0, 1)` by **floor** division of `p/q`, pulling out a `(-1)^a = ±1` sign
  that merges into the coefficient (the reference canonical form). This covers
  negative exponents and `|p| ≥ q` alike, e.g. `(-1)^(-1/5) → -(-1)^(4/5)`,
  `(-1)^(-2/3) → -(-1)^(1/3)`, `(-1)^(-7/5) → (-1)^(3/5)`, `(-1)^(5/4) →
  -(-1)^(1/4)`, `(-8)^(-1/3) → -(-1)^(2/3)/2`. Positive bases keep
  truncation toward zero (residual exponent in `(-1, 1)`), unchanged.
  Negative bases other than `-1` with even `q ≥ 4` (e.g. `(-16)^(1/4)`) are
  still left unevaluated.
- Distributes power over product if the exponent is an integer.
- **Nested rational powers compose for any base when `|inner exponent| < 1`.**
  `(B^r)^s → B^(r·s)` holds on the principal branch for *any* complex `B` when
  `r = p/q` is a non-integer rational with `|p| < q` (then `r·Arg(B)` stays in
  `(-π, π]`, so no branch cut is crossed). This works without a positivity
  assumption: `Sqrt[a^(2/3)] → a^(1/3)`, `Sqrt[Sqrt[a]] → a^(1/4)`. Inner
  exponents with `|p| ≥ q` (e.g. `a^(3/2)`) and symbolic inner exponents (e.g.
  `2^a`) still stay unevaluated.
- **Positive numeric coefficient splits out of a mixed `Times` base** under a
  rational power: `(c·w)^(p/q) → c^(p/q) · w^(p/q)` when `c > 0` is a numeric
  rational/integer/real that fully reduces under `q` and `w` is the symbolic
  residual (valid for any `w` since `Arg(c) = 0`). Combined with the nested-power
  rule this gives `Sqrt[(1/27/a)^(2/3)] → 1/3 (1/a)^(1/3)`. Gated to fire only
  when the coefficient genuinely reduces, so `Sqrt[2 Pi]`, `(4 Pi)^(2/3)`,
  `Sqrt[2 Sqrt[3]]` stay nested.
- **Exact rational base, symbolic exponent** (v0.6146). The reference's normal forms, with `e` any non-numeric exponent (`Pi` and `Sqrt[2]` count) and `c = p/q`, `q > 1`:
  - `(1/n)^e` is `n^(-e)`: `(1/4)^(2 - 2 e)` is `4^(-2 + 2 e)`.
  - `(a/b)^(-r)` is `(b/a)^r` for `a < b`: `(2/3)^(-x)` is `(3/2)^x`, while `(4/3)^(-x)` stays.
  - `b^(c r)` takes an exact `q`-th root `R` of `b`: `4^(-D/2)` is `2^(-D)`, `8^(2 x/3)` is `2^(2 x)`, `(2^70)^(x/2)` is `34359738368^x`.
  - Otherwise, a base whose primes all share one power `K` is `R^(K c r)`: `64^(x/4)` is `2^(3 x/2)`, `1000^(x/2)` is `10^(3 x/2)`. Mixed powers split per prime, with `p^m` kept whole for an integer `m`: `12^(D/2)` is `2^D 3^(D/2)`, `72^(x/2)` is `2^(3 x/2) 3^x`. A root `R < 1` flips to `(1/R)^(-...)`.
  - A negative base splits off `(-1)^(c r)`, which is `I^(p r)` for `q = 2`: `(-1)^(x/2)` is `I^x`, `(-3)^(x/2)` is `I^x 3^(x/2)`. For `q > 2` only when the magnitude rewrites: `(-8)^(x/3)` is `(-1)^(x/3) 2^x`, `(-6)^(x/3)` stays.
  - A squarefree base (`6^(x/2)`, `(2/3)^(x/2)`), an integer coefficient (`4^(2 x)`) and a `Plus` exponent other than under `1/n` (`4^(1 + x/2)`) stay. A bignum base that is not an exact root is left alone where the reference factors it.
  - FeynCalc's Feynman-parameter prefactors meet in this form: `2^(D/2) 4^(-D/2)` is `2^(-D/2)`.
- For `Power[Integer, Rational]` with positive base and positive
  `p/q` exponent, splits the base's prime factorisation into a
  product of distinct-prime powers grouped by reduced effective
  exponent (the reference canonical form). Triggers only when the
  resulting form is strictly more informative -- uniform-exponent
  inputs like `6^(1/3)` and `30^(1/3)` keep the compact form.

```mathematica
In[1]:= Sqrt[45]
Out[1]= 3*Sqrt[5]

In[2]:= (a * b)^2
Out[2]= a^2 * b^2

In[3]:= (-1)^(3/2)
Out[3]= -I

In[4]:= (-1)^(7/4)
Out[4]= -(-1)^(3/4)

In[4b]:= (-1)^(-1/5)        (* floor reduces exponent into [0,1) *)
Out[4b]= -(-1)^(4/5)

In[5]:= 18^(1/3)
Out[5]= 2^(1/3) 3^(2/3)

In[6]:= 12^(1/3)
Out[6]= 2^(2/3) 3^(1/3)

In[7]:= 60^(1/3)            (* 3 and 5 share eff 1/3 -> grouped *)
Out[7]= 2^(2/3) 15^(1/3)

In[8]:= 6^(1/3)             (* uniform exps -> stays *)
Out[8]= 6^(1/3)
```

## Sqrt

Square root.
- `Sqrt[z]`: Internally represented as `Power[z, 1/2]`.

**Integer radicands** — `Power[n, p/q]` for exact integer `n` splits `n` into a
rational coefficient and a `q`-th-power-free residue (`Sqrt[72]` → `6 Sqrt[2]`,
`54^(1/3)` → `3 2^(1/3)`), which needs `n`'s prime factorisation. Small
radicands are trial-divided directly; a residue past that limit goes to the
same factoriser `FactorInteger` uses (Pollard rho / SQUFOF / ECM). Before
2026-07-29 the trial division ran all the way to `sqrt(n)` — up to 3·10⁹ steps
for an `int64` — so a single 18–19 digit non-square radicand cost tens of
seconds (`Sqrt[3141592653589793238]` took 24 s while
`FactorInteger` of the same number took 3.6 ms).

## Mod, Quotient, QuotientRemainder

- `Mod[n, m]`: Remainder of `n/m`.
  - Real operands use the floored quotient: `n - m Floor[n/m]`. `Mod[n, m, d]` uses `n - m Floor[(n - d)/m]`.
  - Since v0.471 Mod also follows the reference beyond integers, rationals and reals:
    - **Complex or Gaussian operands** use the nearest quotient, `n - m Round[(n - d)/m]`, with ties to even: `Mod[3 + 4 I, 2]` → `-1`, `Mod[7, 3 + I, 1/2]` → `I`. Gaussian integers take a direct `Z[i]` route (`nt_gaussian_mod`), in int64 when the modulus components are below 2¹⁶.
      - **Printing difference:** an inexact complex result whose imaginary part rounds to `0.` prints as a Real. Mathilda's arithmetic drops `0.` imaginary parts everywhere; the reference keeps `-0.5 + 0.*I`.
    - **Other exact numerics** use `Floor`: `Mod[Pi, 1]` → `Pi - 3`.
    - **A zero modulus** gives `Indeterminate` with `Mod::indet` (it used to stay unevaluated), on every surface: scalar, packed, and `Compile`, which hands the zero back to the interpreter.
    - **An infinite dividend** gives `Indeterminate` with `Infinity::indet`. A real infinite modulus returns `n` (`Mod[-3, Infinity]` → `-3`), and `ComplexInfinity` gives `Indeterminate`.
    - **Arity**: one argument draws `Mod::argtu`, other counts `Mod::argt`.
  - Packed `Mod[buffer, k]` with a constant `k >= 2` uses a reciprocal-multiply arm (`NDBinaryKernel.whole_i`) with an exact integer correction, so the loop has no hardware division. On a pre-built `Range[10^7]` it went from 0.116 s to 0.037 s (the reference 0.047 s), and the results are bit-identical to the `%` path.
- `Quotient[n, m]`: for **real** `n, m`, `Floor[n/m]` — floored, not truncated or
  rounded, so `Quotient[-5.5, 3.]` is `-2`. The three-argument form
  `Quotient[n, m, d]` is `Floor[(n - d)/m]`. For **complex** `n` or `m` it is
  Gaussian-integer division: the ratio's real and imaginary parts are each
  **rounded to the nearest integer** (ties to even, like `Round[]`), giving the
  quotient of least remainder norm. So `Quotient[10 + 2 I, 3 + I]` is `3`,
  `Quotient[17.5 + 6 I, 1 + 2 I]` is `6 - 6 I`, and `Quotient[5 + 3 I, 2]` is
  `2 + 2 I`. The two branches deliberately differ — an imaginary part changes
  which integer is nearest — and both match the reference.
- `QuotientRemainder[n, m]`: Returns `{Quotient[n, m], Mod[n, m]}`.

Over a list, both run on the machine buffer and keep the result packed. Their
exactness rules differ, and the difference is the point: `Quotient` **narrows** —
its answer is an exact `Integer` whatever it is given, so
`Quotient[{1., 2., 3.}, 3]` is `{0, 0, 1}` — while `Mod` takes the *argument's*
exactness, so `Mod[{1., 2., 3.}, 3]` is `{1., 2., 0.}` and `Mod[{1, 2, 3}, 3]`
is `{1, 2, 0}`. Both match the reference.

Before 2026-08-01 neither could hold an integer buffer, which made them the
break in the chain for integer work: `Mod[Range[10^6] 7919, 1000]` handed 10⁶
separate `Integer`s to whatever came next, and `Union`, `Tally` and
`DeleteDuplicates` downstream ran 176×, 36× and 37× slower than on the same
values packed. `Quotient` was worse than slow — it answered `{0., 1., 1., 2.}`
on a list long enough to pack and the exact `{0, 1, 1, 2}` on one that was not,
the same expression giving different element heads either side of the
250-element threshold.

## IntegerDigits

- `IntegerDigits[n]`: List of the decimal digits of the integer `n`, most
  significant first. The sign of `n` is discarded.
- `IntegerDigits[n, b]`: Base-`b` digits of `n` (requires `b >= 2`).
- `IntegerDigits[n, b, len]`: Pads on the left with leading zeros to length
  `len`. If `n` has more than `len` base-`b` digits, the `len` least-
  significant digits are returned.

**Features**:
- `Protected`, `Listable`. Threading distributes element-wise over a list
  of integers in any argument position, e.g. `IntegerDigits[{6, 7, 2}, 2]`
  and `IntegerDigits[7, {2, 3, 4}]`.
- `IntegerDigits[0]` returns `{0}` (single zero digit). With an explicit
  length, `IntegerDigits[0, b, len]` returns a list of `len` zeros.
- Bases > 10 are allowed; digit values range over `{0, ..., b-1}`.
- Machine integers take a word-sized loop. Bignums in bases up to 62 go
  through GMP's subquadratic `mpz_get_str`, and the digits are written
  straight into a packed `int64` list, which is how the reference returns them:
  `IntegerDigits[3^200000]` takes 0.0035 s (the reference: 0.0038 s; the
  old repeated division took 0.78 s). Larger bases fall back to repeated
  division.
- Messages, as in the reference: `IntegerDigits::ibase` (a base that is not an
  integer greater than 1) and `IntegerDigits::intnm` (a negative length).

```mathematica
In[1]:= IntegerDigits[58127]
Out[1]= {5, 8, 1, 2, 7}

In[2]:= IntegerDigits[58127, 16]
Out[2]= {14, 3, 0, 15}

In[3]:= IntegerDigits[Range[0, 7], 2, 3]
Out[3]= {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}}

In[4]:= IntegerDigits[6345354, 10, 4]
Out[4]= {5, 3, 5, 4}
```

> **Packed arrays.** The one-argument form reads an `int64` buffer directly
> rather than materialising it. The result is ragged — a list per element, of
> differing lengths — so it is an ordinary list either way. The digit nodes and
> list heads are shared (reference-counted), so `IntegerDigits[Range[10^6]]`
> takes 0.36 s (the reference: 0.47 s; before: 1.23 s).

## IntegerLength

- `IntegerLength[n]`: Number of decimal digits in the integer `n`. The sign
  of `n` is discarded.
- `IntegerLength[n, b]`: Number of base-`b` digits in `n` (requires `b >= 2`).

**Features**:
- `Protected`, `Listable`. Threads element-wise over a list of integers in
  either argument position, e.g. `IntegerLength[{1, 10, 100}]` or
  `IntegerLength[8, {2, 3, 4}]`.
- `IntegerLength[0]` is `0` in any base (zero has no significant digits).
- `IntegerLength[n, b]` is effectively an efficient version of
  `Floor[Log[b, |n|]] + 1` -- it never converts through floating-point and
  is exact for arbitrarily large `n`.
- Works for both machine integers and bignums. The fast path uses GMP's
  `mpz_sizeinbase` for bases `2..62` (with an exact verification step for
  non-power-of-2 bases); arbitrary-precision bases fall back to repeated
  division.

```mathematica
In[1]:= IntegerLength[123456789]
Out[1]= 9

In[2]:= IntegerLength[100!, 2]
Out[2]= 525

In[3]:= Table[IntegerLength[100!, n], {n, 2, 20}]
Out[3]= {525, 332, 263, 227, 204, 187, 175, 166, 158, 152, 147, 142, 138, 135, 132, 129, 126, 124, 122}
```

**Arity diagnostics** (`IntegerLength::argt`). Wrong-arity calls emit the
diagnostic and echo the call back unevaluated, matching the reference:

```mathematica
In[4]:= IntegerLength[]
        IntegerLength::argt: IntegerLength called with 0 arguments; 1 or 2 arguments are expected.
Out[4]= IntegerLength[]

In[5]:= IntegerLength[1, 2, 3, 4]
        IntegerLength::argt: IntegerLength called with 4 arguments; 1 or 2 arguments are expected.
Out[5]= IntegerLength[1, 2, 3, 4]
```

**Non-integer-argument diagnostic** (`IntegerLength::int`). A concrete
non-integer numeric (Real, Rational, Complex) at position 1 or 2 emits the
diagnostic and echoes the call back unevaluated; pure symbolic arguments
flow through silently:

```mathematica
In[6]:= IntegerLength[1.1234]
        IntegerLength::int: Integer expected at position 1 in IntegerLength[1.1234].
Out[6]= IntegerLength[1.1234]
```

> **Packed arrays.** The one-argument form runs on an `int64` buffer and
> answers with a packed result.

## IntegerExponent

- `IntegerExponent[n]`: Highest power of `10` that divides `n`; equivalently,
  the number of trailing zeros in the decimal digits of `n`. Equivalent to
  `IntegerExponent[n, 10]`.
- `IntegerExponent[n, b]`: Highest power of `b` that divides `n` (requires
  `b >= 2`). Equivalently, the number of trailing zeros in the base-`b`
  digits of `n`.

**Features**:
- `Protected`, `Listable`. Threads element-wise over a list of integers in
  either argument position, e.g. `IntegerExponent[{10, 100, 1000}]` or
  `IntegerExponent[24, {2, 3, 4, 6}]`.
- Sign of `n` is discarded.
- `IntegerExponent[0, b]` is `Infinity` for any base `b` (every power of `b`
  divides 0).
- Works for both machine integers and bignums. Base 2 uses GMP's
  `mpz_scan1` (lowest-set-bit lookup, O(log n / wordsize)); other bases use
  `mpz_remove` (a single library call that strips all factors of `b`).
- A non-integer numeric base emits `IntegerExponent::ibase`, as in the reference.
- **Machine surfaces.** A packed `int64` list (with a scalar base) runs on a
  buffer kernel and stays packed: 0.021 s for 10^6 elements, against 0.087 s
  in the reference (0.45 s before). A zero element gives `Infinity`, so
  that result is an ordinary list. `IntegerExponent[v, b]` compiles at scalar
  and rank-1 shapes.
- `IntegerExponent[array]` (base 10) maps an int64 packed or visible array
  through the same exact kernel as `IntegerExponent[array, b]`; a 0 element or a
  real array takes the List path (`IntegerExponent[0]` is `Infinity`).

```mathematica
In[1]:= IntegerExponent[1230000]
Out[1]= 4

In[2]:= IntegerExponent[2^10 + 2^7, 2]
Out[2]= 7

In[3]:= IntegerExponent[144, 2]
Out[3]= 4

In[4]:= IntegerExponent[100!, 2]
Out[4]= 97

In[5]:= IntegerExponent[0]
Out[5]= Infinity
```

**Arity diagnostics** (`IntegerExponent::argt`). Wrong-arity calls emit the
diagnostic and echo the call back unevaluated, matching the reference:

```mathematica
In[6]:= IntegerExponent[]
        IntegerExponent::argt: IntegerExponent called with 0 arguments; 1 or 2 arguments are expected.
Out[6]= IntegerExponent[]

In[7]:= IntegerExponent[1, 2, 3, 4]
        IntegerExponent::argt: IntegerExponent called with 4 arguments; 1 or 2 arguments are expected.
Out[7]= IntegerExponent[1, 2, 3, 4]
```

**Non-integer-argument diagnostic** (`IntegerExponent::int`). A concrete
non-integer numeric (Real, Rational, Complex) at position 1 or 2 emits the
diagnostic and echoes the call back unevaluated; pure symbolic arguments
flow through silently:

```mathematica
In[8]:= IntegerExponent[1.123]
        IntegerExponent::int: Integer expected at position 1 in IntegerExponent[1.123].
Out[8]= IntegerExponent[1.123]
```

## DigitSum

- `DigitSum[n]`: Sum of the decimal digits of the integer `n`. Sign of `n` is
  discarded; `DigitSum[0]` returns `0`.
- `DigitSum[n, b]`: Sum of the base-`b` digits of `n` (requires `b >= 2`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

**Examples:**

```mathematica
In[1]:= DigitSum[1234]
Out[1]= 10

In[2]:= DigitSum[255, 16]
Out[2]= 30

In[3]:= DigitSum[{1234, 0, 99}]
Out[3]= {10, 0, 18}
```

**Notes:**
- `DigitSum[n]` is equivalent to `Total[IntegerDigits[n]]`.
- Bignum inputs are supported via GMP.

---

## DigitCount

- `DigitCount[n]`: List of counts of digits `1, 2, ..., 9, 0` in the base-10
  representation of `n`. Sign of `n` is discarded.
- `DigitCount[n, b]`: List of counts of digits `1, 2, ..., b-1, 0` in the
  base-`b` representation of `n` (requires `b >= 2`). Note the ordering:
  digit `0` is **last**, matching the reference's convention.
- `DigitCount[n, b, d]`: Scalar count of digit `d` (an integer with
  `0 <= d < b`) in the base-`b` representation of `n`.

**Features**:
- `Protected`, not `Listable` (as in the reference), but a list in the **first**
  argument threads explicitly: `DigitCount[{1, 0, 1}, 2]` is
  `{{1, 0}, {0, 0}, {1, 0}}`. A non-numeric element emits `DigitCount::int`
  for the whole call, and a list as the base emits `DigitCount::base`.
- `DigitCount[n, b, {d1, d2, ...}]` is the list of the individual counts.
- `0` counts as the single digit 0 (`DigitCount[0]` is
  `{0, 0, 0, 0, 0, 0, 0, 0, 0, 1}`), except in base 2, where it is `{0, 0}`.
  Both match the reference.
- A real `n` counts the digits `RealDigits[n, b]` gives:
  `DigitCount[2.5, 2]` is `{2, 51}`.
- A digit that is not an integer in `[0, b)` emits `DigitCount::digp`.
- Defining property: `Total[DigitCount[n, b]] == IntegerLength[n, b]` for
  any `n != 0`, and `DigitCount[n, b, d] == Count[IntegerDigits[n, b], d]`.
- The 3-argument form supports arbitrary-precision base `b` and digit `d`
  (bignum); the 1/2-argument form caps `b` at `2^20` for the list-returning
  path -- beyond that, `DigitCount::ovfl` fires and the call is left
  unevaluated. Use the 3-argument form for very large bases.

```mathematica
In[1]:= DigitCount[2147, 2, 1]
Out[1]= 5

In[2]:= DigitCount[2147, 2]
Out[2]= {5, 7}

In[3]:= DigitCount[100!]
Out[3]= {15, 19, 10, 10, 14, 19, 7, 14, 20, 30}
```

## FromDigits

- `FromDigits[list]`: Reconstructs an integer from a list of decimal digits,
  most-significant first.
- `FromDigits[list, b]`: Same but in base `b`.
- `FromDigits["string"]`: Constructs an integer from a string of digits.
  Characters `0`-`9`, `a`-`z`, `A`-`Z` represent digit values `0`-`35`.
- `FromDigits["string", b]`: String form in base `b`.
- `FromDigits[{digits, e}]`, `FromDigits[{digits, e}, b]`: the mantissa form
  `RealDigits` returns, an exact rational. A trailing sublist is a recurring
  block: `FromDigits[{{1, 2}, -1}]` is `3/250`, `FromDigits[{{1, 2, {3, 4}}, 2}]`
  is `1222/99`. A list of exponents threads: `FromDigits[{{1, 2}, {3, 4}}]` is
  `{120, 1200}`.
- `FromDigits["MCMXC", "Roman"]`: a Roman numeral (case-insensitive, additive
  forms such as `"IIII"` accepted); an invalid one emits
  `FromRomanNumeral::nrom` and gives `$Failed`.

**Features**:
- `Protected` (intentionally not `Listable`: the first argument *is* a list).
- Inverse of `IntegerDigits` / `IntegerString`. Since `IntegerDigits` discards
  the sign of `n`, `FromDigits[IntegerDigits[n]]` is `Abs[n]`, not `n`.
- Digits in the list (and characters in the string) need *not* be smaller
  than the base; they are carried via Horner's evaluation, matching
  the reference (e.g. `FromDigits[{7, 11, 0, 0, 0, 122}] == 810122`,
  `FromDigits["1A3C"] == 2042`).
- Any integer base is accepted for a list, as in the reference:
  `FromDigits[{1, 2, 3}, 1]` is `6`, `FromDigits[{1, 2, 3}, -10]` is `83`.
- The all-integer case is exact in GMP. Digit lists in bases up to 62 and
  strings in bases up to 36 go through `mpz_set_str` (subquadratic). A packed
  `int64` digit list is read in place (FromDigits is on `AWARE`): 10^5 digits
  take 0.0018 s, against 0.27 s before and 0.0015 s in the reference.
- With a numeric base and symbolic digits the result is the Horner form
  (`FromDigits[{a, b, c}]` is `c + 10 (b + 10 a)`, as in the reference); a
  symbolic base gives the expanded sum `d[n-1] + d[n-2] x + ...`. For a string,
  the base must be an integer; a symbolic base over a string stays unevaluated
  (silently).
- `FromDigits[]` emits `FromDigits::argt`.
- Edge cases: `FromDigits[{}] == 0`, `FromDigits[""] == 0`,
  `FromDigits[{0, 0, 1, 2, 3}] == 123` (leading zeros are inert).

```mathematica
In[1]:= FromDigits[{5, 1, 2, 8}]
Out[1]= 5128

In[2]:= FromDigits[{1, 0, 1, 1, 0, 1, 1}, 2]
Out[2]= 91

In[3]:= FromDigits["1A3C"]
Out[3]= 2042

In[4]:= FromDigits[IntegerDigits[2^100]]
Out[4]= 1267650600228229401496703205376

In[5]:= FromDigits[{a, b, c, d, e}, x]
Out[5]= e + d x + c x^2 + b x^3 + a x^4
```

## IntegerString

- `IntegerString[n]`: String of the decimal digits in the integer `n`. The
  sign of `n` is discarded.
- `IntegerString[n, b]`: Base-`b` digit string (requires `2 <= b <= 36`).
  Digit values `10`-`35` use the lowercase letters `a`-`z`.
- `IntegerString[n, b, len]`: Pads the string on the left with `'0'` to give
  a string of length `len`. If `n` has more than `len` base-`b` digits, the
  `len` least-significant digits are returned (matching the reference).
- `IntegerString[n, "Roman"]` (the name is case-insensitive): `|n|` as a Roman
  numeral -- `IntegerString[2026, "Roman"]` is `"MMXXVI"`, `4000` is `"MMMM"`,
  `0` is `"N"`. A length argument does not pad it. `|n| >= 5000` needs
  overbarred numerals, a typeset form rather than a string, and stays
  unevaluated. Written from the same numeral table `FromDigits[s, "Roman"]`
  reads, so the two round-trip.

**Features**:
- `Protected`, `Listable`. Threads element-wise over a list of integers in
  any argument position, e.g. `IntegerString[Range[0, 7], 2, 3]` returns
  a length-8 list of zero-padded binary strings.
- Inverse of `FromDigits` for the integer case:
  `FromDigits[IntegerString[n, b], b] == n` for any non-negative `n`.
- Defining property: `StringLength[IntegerString[n, b]] == IntegerLength[n, b]`
  for `n != 0`.
- `IntegerString[0]` returns `"0"`; `IntegerString[n, b, 0]` returns `""`.
- Works seamlessly for machine integers and arbitrary-precision bignums --
  the digit rendering goes through GMP's `mpz_get_str` in a single call,
  so even `IntegerString[100!, 36]` is essentially free.

```mathematica
In[1]:= IntegerString[17651, 2]
Out[1]= "100010011110011"

In[2]:= IntegerString[50!, 16]
Out[2]= "49eebc961ed279b02b1ef4f28d19a84f5973a1d2c7800000000000"

In[3]:= IntegerString[50!, 36]
Out[3]= "4q7eyp9zizmtqt0648txt4fm720cc1s00000000000"

In[4]:= IntegerString[Range[0, 7], 2, 3]
Out[4]= {"000", "001", "010", "011", "100", "101", "110", "111"}

In[5]:= IntegerString[12345, 10, 3]
Out[5]= "345"
```

**Arity diagnostics** (`IntegerString::argb`). Wrong-arity calls emit the
diagnostic and echo the call back unevaluated, matching the reference:

```mathematica
In[6]:= IntegerString[]
        IntegerString::argb: IntegerString called with 0 arguments; between 1 and 3 arguments are expected.
Out[6]= IntegerString[]
```

**Non-integer-argument diagnostic** (`IntegerString::int`). A concrete
non-integer numeric (Real, Rational, Complex) at position 1 or 3 emits the
diagnostic and echoes the call back unevaluated; pure symbolic arguments
flow through silently:

```mathematica
In[7]:= IntegerString[11.3423]
        IntegerString::int: Integer expected at position 1 in IntegerString[11.3423].
Out[7]= IntegerString[11.3423]
```

**Invalid-base diagnostic** (`IntegerString::basf`). Bases outside the
`[2, 36]` range (integer or otherwise) trigger the diagnostic and the call
is left unevaluated:

```mathematica
In[8]:= IntegerString[10, 50]
        IntegerString::basf: 50 is not a valid base for IntegerString in IntegerString[10, 50]; the base must be an integer between 2 and 36.
Out[8]= IntegerString[10, 50]
```

## RealDigits

- `RealDigits[x]`: List of the base-10 digits of the approximate real
  number `x`, together with the integer exponent: the result is
  `{digits, exp}` where the first element of `digits` is the coefficient
  of `10^(exp - 1)`.
- `RealDigits[x, b]`: Same in base `b` (requires integer `b >= 2`).
- `RealDigits[x, b, len]`: Exactly `len` digits.
- `RealDigits[x, b, len, n]`: `len` digits, first one the coefficient of
  `b^n`; the returned exponent is `n + 1`.

**Features**:
- `Protected`, `Listable`. Threads over lists in any argument position.
- Works for `Integer`, `BigInt`, `Rational`, machine `Real`, and (under
  `USE_MPFR`) arbitrary-precision `MPFR` numbers. The sign of `x` is
  discarded.
- For integers and rationals with terminating base-`b` expansions, the
  digit list is flat. For rationals with non-terminating expansions, the
  list ends in a nested list giving the recurring block:
  `RealDigits[19/7]` returns `{{2, {7, 1, 4, 2, 8, 5}}, 1}`.
- For inexact (machine or MPFR) reals, the default `len` is
  `Round[Precision[x] / Log10[b]]`. Requesting more digits than the
  precision allows produces `Indeterminate` at the LSB end. The digits
  themselves use the canonical round-to-nearest representation supplied
  by MPFR -- so `RealDigits[123.55555]` returns the literal decimal
  digits, not the binary IEEE tail.
- Symbolic numeric constants such as `Pi`, `E`, and `GoldenRatio` are
  numericalized to MPFR at the matching precision when an explicit `len`
  is given. `RealDigits[Pi]` (no `len`) is left unevaluated.
- `RealDigits[0]` returns `{{0}, 0}`. `RealDigits[0.]` returns
  `{{0}, -Floor[Accuracy[0.]]}` — `{{0}, -323}` for machine precision
  (`Accuracy[0.] ≈ 323.607`), and `{{0}, -p}` for an MPFR zero
  `0``p` of precision `p` digits.
- Bases must be integers `>= 2`. Non-integer (Real / Rational) bases
  trigger `RealDigits::ibase` and leave the call unevaluated.
  Non-integer-base expansions (e.g. `GoldenRatio`) are not yet supported.
- `FromDigits` can be used as the inverse of `RealDigits` for the
  integer / terminating-rational case.

```mathematica
In[1]:= RealDigits[123.55555]
Out[1]= {{1, 2, 3, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0}, 3}

In[2]:= RealDigits[Pi, 10, 25]
Out[2]= {{3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3}, 1}

In[3]:= RealDigits[19/7]
Out[3]= {{2, {7, 1, 4, 2, 8, 5}}, 1}

In[4]:= RealDigits[5.635, 10, 20]
Out[4]= {{5, 6, 3, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Indeterminate, Indeterminate, Indeterminate, Indeterminate}, 1}

In[5]:= RealDigits[Pi, 10, 20, -5]
Out[5]= {{9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3}, -4}

In[6]:= RealDigits[1.234, 2, 15]
Out[6]= {{1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1}, 1}

In[7]:= {RealDigits[1000.], RealDigits[0.999]}
Out[7]= {{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 4}, {{9, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 0}}
```

## MantissaExponent

- `MantissaExponent[x]`: Returns `{m, e}` such that `x == m * 10^e` and
  `1/10 <= |m| < 1` (or `{0, 0}` when `x` is zero).
- `MantissaExponent[x, b]`: Same in base `b`; the mantissa then lies in
  `1/b <= |m| < 1`.

**Features**:
- `Protected`, `Listable`. Threads over lists in any argument position.
- Works for `Integer`, `BigInt`, `Rational`, machine `Real`, and (under
  `USE_MPFR`) arbitrary-precision `MPFR` numbers. The sign of `x` carries
  through to the mantissa.
- For exact inputs the mantissa is an exact `Rational` of the form
  `x / b^e`; for inexact inputs the mantissa is returned with the same
  precision as `x`.
- For a machine real the exponent is decided on the double's exact value, so a
  value just above a power of the base is on the right side of it
  (`MantissaExponent[0.01]` is `{0.1, -1}`, since the double 0.01 is
  0.0100000000000000002...), and the mantissa is rounded once.
- Bases must be integers `>= 2`. Non-integer bases leave the call
  unevaluated (general `Real` / symbolic bases are not yet supported).
- Complex arguments emit `MantissaExponent::realx` and leave the call
  unevaluated. Symbolic (non-numeric) arguments are left unevaluated
  silently.

```mathematica
In[1]:= MantissaExponent[3.4 10^30]
Out[1]= {0.34, 31}

In[2]:= MantissaExponent[456.1414]
Out[2]= {0.456141, 3}

In[3]:= MantissaExponent[123451]
Out[3]= {123451/1000000, 6}

In[4]:= MantissaExponent[1027, 2]
Out[4]= {1027/2048, 11}

In[5]:= MantissaExponent[2^100, 2]
Out[5]= {1/2, 101}

In[6]:= MantissaExponent[N[Pi, 30]]
Out[6]= {0.3141592653589793238462643383278, 1}

In[7]:= MantissaExponent[-3/2]
Out[7]= {-3/20, 1}
```

## RealExponent

- `RealExponent[x]`: Returns `Log[10, |x|]` -- the base-10 real exponent
  of `x`.
- `RealExponent[x, b]`: Returns `Log[b, |x|]` in the specified base `b`.

**Features**:
- `Protected`, `Listable`. Threads over lists in any argument position.
- Accepts `Integer`, `BigInt`, `Rational`, machine `Real`, and (under
  `USE_MPFR`) arbitrary-precision `MPFR` inputs.  Symbolic numeric
  arguments (`Pi`, `E`, `EulerGamma`, `Catalan`, `GoldenRatio`, `Degree`,
  or any numeric-valued composite such as `Pi^Pi` or `1/Pi`) are
  numericalized at the combined working precision before computation.
  Plain symbols with no numeric value are left unevaluated.
- Output is a machine `Real` unless one of the inputs already carries
  MPFR precision, in which case the result is `MPFR` at the higher of
  the input precisions.  This matches the reference's contagion: an
  explicit `N[..., p]` lifts the exponent to the same precision.
- The base must be a real number `> 1`; non-positive, `<= 1`, or complex
  bases emit `RealExponent::ibase` and leave the call unevaluated.
- Complex arguments with non-zero imaginary part emit
  `RealExponent::realx` and leave the call unevaluated.
- Sign of `x` is discarded.
- Zero handling (reference-compatible):
  - Exact zero (Integer 0, BigInt 0, Rational 0/n) → `-Infinity`.
  - Machine `0.` → `Log[b, $MinMachineNumber]` (`≈ -307.65` for base 10).
  - MPFR `0``p` of `p` digits → `-p / Log10[b]` (`-p` for base 10).

```mathematica
In[1]:= RealExponent[123.456]
Out[1]= 2.09151

In[2]:= RealExponent[123.456, 2]
Out[2]= 6.94785

In[3]:= RealExponent[N[Pi, 32]]
Out[3]= 0.497149872694133854351268288290899

In[4]:= RealExponent[Pi, E]
Out[4]= 1.14473

In[5]:= RealExponent[987654321/123456789]
Out[5]= 0.90309

In[6]:= RealExponent[{1, 2, 3, 4, 5}]
Out[6]= {0.0, 0.30103, 0.477121, 0.60206, 0.69897}

In[7]:= Table[RealExponent[Pi, b], {b, {2, 3, 5, 7, 10}}]
Out[7]= {1.6515, 1.04198, 0.711261, 0.588275, 0.49715}

In[8]:= RealExponent[0]
Out[8]= -Infinity
```

## Factorial (!)

Gives the factorial of an integer or half-integer.
- `n!` or `Factorial[n]`

**Features**:
- `Protected`, `Listable`, `NumericFunction`.
- Evaluates exactly for positive integers up to $20!$.
- Yields `ComplexInfinity` for negative integers.
- Supports half-integers utilizing factors of $\sqrt{\pi}$ recursively.
- Supports trailing `!` parsed natively as a postfix operator.

```mathematica
In[1]:= 5!
Out[1]= 120

In[2]:= (1/2)!
Out[2]= Sqrt[Pi]/2

In[3]:= Factorial[0]
Out[3]= 1
```

## FactorialPower

Falling factorial. `FactorialPower[n, k]` = $n (n-1) (n-2) \cdots (n - k + 1)$.
- `FactorialPower[n, k]`

**Features**:
- `Protected`, `Listable`, `NumericFunction`.
- For non-negative integer $n$ and non-negative integer $k$: exact GMP product.
- For symbolic $n$ with concrete $k \le 32$: expands to an explicit product
  `Times[n, n - 1, ..., n - k + 1]` so `Expand` and `D` can act on it.
- Equivalent to $n! / (n - k)!$ when both arguments are non-negative integers.
- For a **non-integer numeric** $k$ (under `N`), evaluates the generalised form
  $\Gamma(n+1)/\Gamma(n-k+1)$ (real, complex, arbitrary precision); exact
  non-integer $k$ stays symbolic.
- Used as the symbolic-order derivative of a power: $D[x^n, \{x, k\}] = \mathrm{FactorialPower}[n, k]\, x^{n-k}$.

```mathematica
In[1]:= FactorialPower[5, 3]
Out[1]= 60

In[2]:= FactorialPower[n, 3]
Out[2]= n (-2 + n) (-1 + n)

In[3]:= FactorialPower[n, 0]
Out[3]= 1

In[4]:= D[x^n, {x, k}]
Out[4]= x^(-k + n) FactorialPower[n, k]
```

- **v0.477**: the step form `FactorialPower[x, n, h]` = x (x - h) … (x - (n - 1) h)
  for a positive integer `n`, 1/((x + h) … (x + |n| h)) for a negative one, and
  h^n Gamma[x/h + 1]/Gamma[x/h - n + 1] otherwise, evaluated when `x` and `h` are
  numeric (`FactorialPower[5, 3, 2] = 15`, `FactorialPower[5, 1/2, 2] = 15
  Sqrt[Pi/2]/8`, `FactorialPower[-2.5, 1.5, 0.5] = ComplexInfinity`); a symbolic
  `x` or `h` stays. A negative integer `k` with numeric `n` now evaluates too
  (`FactorialPower[5, -2] = 1/42`).

## Binomial

Gives the binomial coefficient $\binom{n}{m}$, generalised via
$\Gamma(n+1)/(\Gamma(m+1)\,\Gamma(n-m+1))$.

- `Binomial[n, m]`

**Features**:
- `Protected`, `Listable`, `NumericFunction`.
- Exact integer/integer path uses GMP (`mpz_bin_uiui` when `n` is a machine
  word -- ~1 us for `Binomial[355, 150]`, against ~7 us for the general
  `mpz_bin_ui` it used before; `mpz_bin_ui` for larger `n`), including the
  Pascal extension `Binomial[n, m] = (-1)^m Binomial[m-n-1, m]` for
  negative `n`. When both are negative the Gamma quotient's poles cancel for
  `m ≤ n`: `Binomial[n, m] = (-1)^(n-m) Binomial[-m-1, n-m]`, so
  `Binomial[-1, -1] = 1` and `Binomial[-2, -5] = -4`; `m > n` gives `0`.
  The smaller of `m` and `n - m` is used, so `Binomial[2^70 + 1, 2^70]`
  evaluates.
- Machine-precision real branch via `tgamma`. Fires when either operand
  is a machine real; the other may be any real-numeric leaf (integer,
  bigint, rational, MPFR) — so `Binomial[2.5, 1/5] → 1.34885`. If an
  operand is symbolic the branch declines and the form stays symbolic
  (`Binomial[2.5, x]`) rather than reading the sibling as `0`.
- Symmetric identity: when `n - m` simplifies to a non-negative
  integer `k ≤ 32` (`k ≤ 5` when `n` is not an explicit number, as in the
  reference: `Binomial[a + 5, a]` expands, `Binomial[a + 6, a]` stays), reduces to `Binomial[n, k]` and expands as a
  falling-factorial polynomial. This catches `Binomial[n, n - 1] → n`,
  `Binomial[9/2, 7/2] → 9/2`, `Binomial[n + 1, n - 1] → n (n + 1)/2`,
  etc.
- Concrete non-negative integer `m ≤ 32` (`m ≤ 5` for a symbolic `n`) with
  any other `n` (symbolic, rational, complex, …) expands to the falling-factorial polynomial
  `n (n-1) (n-2) ... (n-m+1) / m!`, which the `Times`/`Plus` folders
  then simplify.
- Exact rationals off the integer lattice. If `m` or `n - m` is a negative
  integer and `n` is not, the result is `0` (`Binomial[1/3, 4/3] = 0`). If
  `n`, `m` and `n - m` are all multiples of `1/2`, every Gamma is a rational
  multiple of a power of `Sqrt[Pi]` and the quotient closes:
  `Binomial[0, 1/2] = 2/Pi` and `Binomial[5, 1/2] = 512/(63 Pi)`. When `n` is a
  negative integer the result is `ComplexInfinity` (`Binomial[-2, 1/2]`). Any
  other rational pair stays symbolic (`Binomial[1/3, 1/2]`).
- Arbitrary-precision real branch via MPFR (`mpfr_gamma`), evaluating the
  Gamma quotient at the working precision. It leaves the exact rational
  form symbolic — `Binomial[7/3, 1/5]` stays `Binomial[7/3, 1/5]`, as in
  the reference — but supplies a value under `N[Binomial[7/3, 1/5], p]`.
  Placed after the integer-`m` branch so an integer `m` keeps the exact
  falling-factorial form. The working precision follows the usual
  contagion rule (minimum precision among the inexact operands, floored
  at machine 53).
- Complex numeric branch: fires when the computation is numeric — an
  inexact operand is present (a machine or MPFR real, or a `Complex[..]`
  with an inexact part) — **and** a complex operand is present, so the
  real-only machine branch has declined. The Gamma quotient is built as an
  expression and evaluated, reusing `Gamma`'s complex kernels (machine
  Lanczos, and the arbitrary-precision Spouge path under MPFR) and the
  complex-arithmetic folders; any exact Gaussian sibling (`1 + I`,
  `7 - 3 I`) is carried along by the `Times`/`Plus` numeric contagion as the
  quotient folds. So `Binomial[1 + I, 5.]`, `Binomial[2. + I, 7 - 3 I]`,
  `N[Binomial[1/2 + I/3, 1/4]]` and `N[Binomial[1/2 + I/3, 1/4], 25]` all
  evaluate. The result is accepted only if it is numeric, so a pole or a
  symbolic operand leaves the form symbolic. A pair of **exact** Gaussians
  with no inexact operand (`Binomial[1 + I, 2 + I]`) stays symbolic.

```mathematica
In[1]:= Binomial[10, 3]
Out[1]= 120

In[2]:= Binomial[8.5, -4.2]
Out[2]= 0.0000604992

In[3]:= Binomial[9/2, 7/2]
Out[3]= 9/2

In[4]:= Binomial[n, 4]
Out[4]= 1/24 n (-3 + n) (-2 + n) (-1 + n)

In[5]:= Binomial[n, n - 1]
Out[5]= n

In[6]:= Binomial[1 + I, 5]
Out[6]= -1/12 - I/12

In[7]:= Binomial[0, 1]
Out[7]= 0

In[8]:= Binomial[7/3, 1/5]
Out[8]= Binomial[7/3, 1/5]

In[9]:= N[Binomial[7/3, 1/5]]
Out[9]= 1.33313

In[10]:= N[Binomial[7/3, 1/5], 25]
Out[10]= 1.3331254244650286522359229

In[11]:= Binomial[2.5, 1/5]
Out[11]= 1.34885

In[12]:= Binomial[1/2 + I/3, 1/4]
Out[12]= Binomial[1/2 + I/3, 1/4]

In[13]:= N[Binomial[1/2 + I/3, 1/4]]
Out[13]= 1.08987 + 0.0929283 I

In[14]:= N[Binomial[1/2 + I/3, 1/4], 25]
Out[14]= 1.0898678407199392604353272 + 0.092928304677202434313313055 I

In[15]:= Binomial[1 + I, 5.]
Out[15]= -0.0833333 - 0.0833333 I

In[16]:= Binomial[2. + I, 7 - 3 I]
Out[16]= -75.4683 + 106.815 I
```

## Multinomial

The multinomial coefficient $(n_1 + n_2 + \dots)!/(n_1!\,n_2!\cdots)$.

- `Multinomial[n1, n2, ...]`

**Features**:
- `Protected`, `Listable`, `NumericFunction`, `Orderless`.
- Evaluated as the product of binomials over the partial sums of the sorted
  arguments, `Binomial[n1 + n2, n2] Binomial[n1 + n2 + n3, n3] ...`, so it
  inherits `Binomial`'s closed forms and pole conventions:
  `Multinomial[-5, 2, 1] = -12`, `Multinomial[1/2, -1/2] = 2/Pi`.
- Exact zeros drop out (`Multinomial[a, b, 0]` is `Binomial[a + b, b]`).
  `Multinomial[]` and `Multinomial[x]` are `1`.
- Two arguments always give `Binomial[n1 + n2, n2]`.
- Three or more give the product only when at most one argument is not an
  explicit number and at most one binomial stays unevaluated:
  `Multinomial[2, 3, 4, a] = 1260 Binomial[9 + a, a]`,
  `Multinomial[1, 1/3, 1/5] = 23/15 Binomial[8/15, 1/3]`. Otherwise the call
  stays (`Multinomial[2, a, b]`, `Multinomial[1, E, Pi]`,
  `Multinomial[1/5, 1/3, 1/2]`). These rules were measured against the
  reference.
- `Compile[]` lowers it to the same product of binomials, at scalar and array
  shapes. A compiled body sums in the written order, while the interpreter sums
  in the sorted order, so a machine-real answer can differ in the last bit.

```mathematica
In[1]:= Multinomial[2, 3, 4, 5]
Out[1]= 2522520

In[2]:= Multinomial[a, 2]
Out[2]= 1/2 (1 + a) (2 + a)

In[3]:= Multinomial[1.5, 2, a]
Out[3]= 4.375 Binomial[3.5 + a, a]
```

## Fibonacci

- `Fibonacci[n]`: the `n`th Fibonacci number `F_n`.
- `Fibonacci[n, x]`: the `n`th Fibonacci polynomial `F_n(x)`, the coefficient
  of `t^n` in `t / (1 - x t - t^2)`.
- Attributes: `Listable`, `NumericFunction`, `Protected`.

Evaluation:
- Exact integer `n`: `Fibonacci[n]` uses GMP **fast doubling** (`O(log n)`);
  `Fibonacci[n, x]` builds the polynomial from the recurrence
  `F_k = x F_{k-1} + F_{k-2}` (Expand-ed each step). Negative orders use
  `F_{-n} = (-1)^{n+1} F_n` (and likewise for the polynomials).
- Inexact / complex order (`Real`, MPFR, or `Complex`): the generalized closed
  forms — `Fibonacci[n] = (φ^n - Cos[π n] φ^{-n}) / Sqrt[5]` with
  `φ = GoldenRatio`, and `Fibonacci[n, x] = (β^n - Cos[π n] β^{-n}) / Sqrt[x²+4]`
  with `β = (x + Sqrt[x²+4]) / 2` — are numericalized at the precision carried
  by the inputs (machine or MPFR).
- Exact non-integer order with a numeric `x` (`Fibonacci[n, x]` only): the same
  closed form is evaluated *exactly*. An inexact `x` numericalizes it
  (`Fibonacci[1/2, 3.2] = 0.494833`); an exact `x` keeps the result only when it
  collapses to a number (`Fibonacci[1/2, 0] = 1/2`, since `F_n(0) = (1 - Cos[π n])/2`),
  otherwise the call stays unevaluated.
- Symbolic order, symbolic `x`, or a one-argument exact non-integer order with
  no `N` applied stays unevaluated.

Derivatives (native, via `src/calculus/deriv.c`):
- `D[Fibonacci[n], n]` and both partials of `D[Fibonacci[n, x], ...]` are
  provided, e.g.
  `D[Fibonacci[n, x], x] = (2 n Fibonacci[n-1, x] + (n-1) x Fibonacci[n, x]) / (4 + x²)`.

```
In[1]:= Table[Fibonacci[n], {n, 10}]
Out[1]= {1, 1, 2, 3, 5, 8, 13, 21, 34, 55}
In[2]:= Fibonacci[7, x]
Out[2]= 1 + 6 x^2 + 5 x^4 + x^6
In[3]:= Fibonacci[5.8, 3]
Out[3]= 283.483
In[4]:= N[Fibonacci[15/17], 50]
Out[4]= 0.95651991392431122508582263427692298648606969012061
In[5]:= Fibonacci[1/2, 0]
Out[5]= 1/2
In[6]:= Fibonacci[1/2, 3.2]
Out[6]= 0.494833
```

## LucasL

- `LucasL[n]`: the `n`th Lucas number `L_n`.
- `LucasL[n, x]`: the `n`th Lucas polynomial `L_n(x)`, the coefficient of
  `t^n` in `(2 - t x) / (1 - x t - t^2)`.
- Attributes: `Listable`, `NumericFunction`, `Protected`.

Evaluation:
- Exact integer `n`: `LucasL[n]` uses GMP **fast doubling** (`O(log n)`) of the
  Fibonacci pair, `L_m = 2 F_{m+1} - F_m`; `LucasL[n, x]` builds the polynomial
  from the recurrence `L_k = x L_{k-1} + L_{k-2}` with `L_0 = 2`, `L_1 = x`
  (Expand-ed each step). Negative orders use `L_{-n} = (-1)^n L_n` (and likewise
  for the polynomials).
- Inexact / complex order (`Real`, MPFR, or `Complex`): the generalized closed
  forms — `LucasL[n] = φ^n + Cos[π n] φ^{-n}` with `φ = GoldenRatio`, and
  `LucasL[n, x] = β^n + Cos[π n] β^{-n}` with `β = (x + Sqrt[x²+4]) / 2` — are
  numericalized at the precision carried by the inputs (machine or MPFR).
- Symbolic order (or exact non-integer order with no `N` applied) stays
  unevaluated.

Derivatives (native, via `src/calculus/deriv.c`):
- `D[LucasL[n], n]` and both partials of `D[LucasL[n, x], ...]` are provided,
  e.g. `D[LucasL[n, x], x] = (n (2 LucasL[n-1, x] + x LucasL[n, x])) / (4 + x²)`.

```
In[1]:= Table[LucasL[n], {n, 10}]
Out[1]= {1, 3, 4, 7, 11, 18, 29, 47, 76, 123}
In[2]:= LucasL[7, x]
Out[2]= 7 x + 14 x^3 + 7 x^5 + x^7
In[3]:= LucasL[-11.]
Out[3]= -199.
In[4]:= N[LucasL[11/3], 50]
Out[4]= 5.9239626529619554101356978621940126287519855422362
```

## Rational

Represents a rational number.
- `Rational[n, d]`

**Features**:
- Automatically simplifies to lowest terms (e.g. `Rational[15, 5]` evaluates to `3`, `Rational[2, 4]` evaluates to `Rational[1, 2]`).
- Returns `Indeterminate` when `n` and `d` are both `0` (e.g. `Rational[0, 0]`).
- Returns `ComplexInfinity` when `n` is non-zero and `d` is `0` (e.g. `Rational[1, 0]`).

```mathematica
In[1]:= Rational[15, 5]
Out[1]= 3
```

## Re, Im, ReIm, Abs, Sign, Conjugate, Arg

Complex number functions.
- `Re[z]`: Real part. Folds `Re[f[z]] -> f[z]` when `f` is one of `Re`,
  `Im`, `Abs`, `Arg` (all real-valued by construction).
- `Im[z]`: Imaginary part. Folds `Im[f[z]] -> 0` for the same set of `f`.
- `ReIm[z]`: Returns `{Re[z], Im[z]}`; inherits the same folding.
- `Abs[z]`: Magnitude. Numeric arguments (`Integer`, `BigInt`, `Real`,
  arbitrary-precision `MPFR`, `Rational`) reduce to a concrete value at the
  input precision; complex arguments use `Sqrt[Re^2 + Im^2]`. When `z` is
  `Complex[MPFR, MPFR]` (or `Complex` with any MPFR component), `Abs`
  takes a direct `mpfr_hypot` fast path at the working precision.
- `Sign[z]`: For real numeric `z`, returns the exact integer `-1`, `0`, or
  `1` (regardless of whether the input is `Integer`, `BigInt`, `Real`,
  `MPFR`, or `Rational`). For nonzero numeric complex `z`, returns
  `z / Abs[z]`. `Complex[MPFR, MPFR]` goes through a direct
  `mpfr_div`-on-components path at MPFR precision. `Listable`.
- `Conjugate[z]`: Complex conjugate.  Folds `Conjugate[Conjugate[z]] -> z`
  (involution) and `Conjugate[f[..]] -> f[..]` when `f` is one of `Re`,
  `Im`, `Abs`, `Arg` (all real-valued by construction).
- `Arg[z]`: Phase angle. Pure-real MPFR returns symbolic `0` or `Pi`;
  `Complex[MPFR, MPFR]` evaluates via `mpfr_atan2` at the working
  precision.  A **root-of-unity-like constant** — a product of real-rational
  powers and the imaginary unit, e.g. `(-1)^(1/6)`, `-(-1)^(5/6)`, `(-8)^(1/3)`,
  `I (-1)^(1/4)` — returns the exact rational multiple of `Pi`
  (`Arg[(-1)^(1/6)] = Pi/6`, `Arg[-(-1)^(1/6)] = -5 Pi/6`), reduced into
  `(-Pi, Pi]`.  The argument adds mod `2 Pi` over the factors (`arg = ([c<0] +
  Σ eᵢ[bᵢ<0]) Pi`), so the generic `a + b I` re/im split that cannot recover it
  is bypassed; a free symbol or a non-root-of-unity Gaussian stays unevaluated.

**Exact numeric constants.** `Re`, `Im`, `ReIm`, `Abs`, `Sign` and `Arg` also
evaluate on closed-form constants such as `Sqrt[2]`, `Pi - 4`, `Sin[10]`,
`E + I Pi` or `1 + I Sqrt[2]`, as the reference does: `Abs[Pi - 4] = 4 - Pi`,
`Sign[Pi - 3] = 1`, `Re[E + I Pi] = E`, `Sign[I Pi] = I`, `Arg[-Sqrt[2]] = Pi`,
`Arg[-1 + I Sqrt[2]] = Pi - ArcTan[Sqrt[2]]`, `Abs[1 + I Sqrt[2]] = Sqrt[3]`.
The real and imaginary parts come from the structural `a + b I` split, and the
**sign** of each part is *certified* by a rigorous FLINT/Arb ball evaluation
(precision ladder 64 → 256 → 1024 bits; multi-valued `Log`/`Power`/inverse
trig are evaluated only off their branch cuts), never by a rounded `N[]`. An
expression that is exactly zero but does not autosimplify — e.g.
`Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]` — has a ball that straddles 0 at
every precision, so the head is left unevaluated rather than guessing a sign.
The supported constant language is: rationals and Gaussian rationals; `Pi`,
`E`, `EulerGamma`, `Catalan`, `GoldenRatio`, `Degree`, `Khinchin`, `Glaisher`;
`Plus`, `Times`, `Power`, `Sqrt`, `Exp`, `Log`, the trig, hyperbolic and
inverse trig/hyperbolic functions, `Gamma`, `Erf`, `Erfc`, `Abs`, `Re`, `Im`,
`Conjugate`. Anything outside it (a free symbol, an inexact leaf, another head)
leaves the call unevaluated, as before.

```mathematica
In[1]:= Re[1/(2 + 3 I)]
Out[1]= 2/13

In[2]:= Im[Log[-1]]
Out[2]= Pi

In[3]:= ReIm[(2 + I)^3]
Out[3]= {2, 11}

In[4]:= Abs[-9223372036854775807 - 1]
Out[4]= 9223372036854775808

In[5]:= Sign[3 + 4 I]
Out[5]= 3/5 + 4/5*I

In[6]:= Conjugate[(2 + I)/(1 - 3 I)]
Out[6]= -1/10 - 7/10*I

In[7]:= Arg[(-1)^(1/6)]
Out[7]= 1/6 Pi

In[8]:= {Abs[Pi - 4], Sign[Sqrt[2]], Im[Pi - 4], Arg[Pi - 4]}
Out[8]= {4 - Pi, 1, 0, Pi}

In[9]:= {Re[E + I Pi], Sign[I Pi], Arg[-1 + I Sqrt[2]]}
Out[9]= {E, I, Pi - ArcTan[Sqrt[2]]}

In[10]:= Abs[Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]]   (* exactly 0: undecidable sign, left alone *)
Out[10]= Abs[Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]]
```

## ComplexExpand

`ComplexExpand[expr]` rewrites `expr` into explicit real and imaginary parts,
assuming every free symbol is **real**.  `ComplexExpand[expr, {x1, x2, ...}]`
instead treats variables matching any of the `xi` (which may be patterns) as
**complex**, breaking each into `Re[xi]` and `Im[xi]`.

- Decomposes recursively through `Plus`, `Times`, `Power` (Cartesian for
  integer exponents, a polar `Abs`/`Arg` master formula otherwise), `Exp`,
  `Log`, the circular and hyperbolic functions and their inverses (rewritten
  through their logarithmic forms), and the `Re`/`Im`/`Abs`/`Arg`/
  `Conjugate`/`Sign`/`ReIm` heads.
- The result is written `re + I im` with each part expanded and `I` kept
  outside the imaginary part, as in the reference (v0.6136):
  `ComplexExpand[(1 + 2 I) x + (3 + 4 I) y]` is `x + 3 y + I (2 x + 4 y)`.
- Option `TargetFunctions` selects the output basis:
  `-> {Re, Im}` (default, Cartesian), `-> {Abs, Arg}` (polar), or
  `-> Conjugate`.
- Threads over `List`, and over equations, inequalities, and logic heads
  (`Equal`, `Less`, `And`, `Or`, ...); verified complex identities such as
  `ComplexExpand[Re[z] == (z + Conjugate[z])/2, z]` collapse to `True`.
- `Protected`.  Calls with zero or more than two positional arguments emit
  `General::argct` and stay unevaluated.

Examples:
```
In[1]:= ComplexExpand[Sin[x + I y]]
Out[1]= Sin[x] Cosh[y] + I Cos[x] Sinh[y]

In[2]:= ComplexExpand[Re[z^2], {z}]
Out[2]= -Im[z]^2 + Re[z]^2

In[3]:= ComplexExpand[Re[z^2], {z}, TargetFunctions -> Conjugate]
Out[3]= 1/2 (z^2 + Conjugate[z]^2)

In[4]:= ComplexExpand[Tan[x + I y]]
Out[4]= Sin[2 x]/(Cos[2 x] + Cosh[2 y]) + I Sinh[2 y]/(Cos[2 x] + Cosh[2 y])
```

---

## Interval arithmetic

`Interval[{min, max}]` represents the closed range of real values between `min`
and `max`. The multi-argument form
`Interval[{min1, max1}, {min2, max2}, ...]` is the union of the ranges. An
`Interval` is an ordinary expression with head `Interval` (like `Complex[re,im]`
and `Rational[n,d]`) and participates in arithmetic and elementary functions —
this is *interval arithmetic*, where every result is a rigorous enclosure of the
true value.

**Representation and canonical form.** Each argument is a two-element `List`
`{lo, hi}`; the canonicalizer orders each pair so `lo ≤ hi`, sorts the pairs by
lower endpoint, and merges overlapping or touching pairs (closed-interval union,
so `{1,2}` and `{2,3}` merge to `{1,3}`). `Interval[x]` for a single number is
the point interval `Interval[{x, x}]`. Endpoints may be exact
(Integer/Rational), `±Infinity`, symbolic numerics (`Sin[2]`, `Pi`), or
approximate (Real/MPFR). Non-comparable symbolic endpoints are left in place.

**Exactness.** Exact endpoints stay exact/symbolic; only approximate endpoints
go numeric. Endpoint operations are carried out by the ordinary (correctly
rounded) kernels, then an inexact result is nudged one ULP *outward* (lower
bounds down, upper bounds up) so the enclosure stays rigorous — no new
dependency is required (the widening uses `nextafter` / `mpfr_nextbelow`
/`mpfr_nextabove`).

**Arithmetic.** `Plus`, `Times`, `Power` (integer and positive-real exponents),
`Divide`, and `Subtract` thread through intervals, including disjoint-union
results:

```
In[1]:= Interval[{1, 6}] + Interval[{0, 2}]
Out[1]= Interval[{1, 8}]

In[2]:= Interval[{-2, 5}]^2
Out[2]= Interval[{0, 25}]

In[3]:= 1/Interval[{-2, 5}]
Out[3]= Interval[{-Infinity, -1/2}, {1/5, Infinity}]

In[4]:= Interval[{-1, 2}] Interval[{3, 4}]
Out[4]= Interval[{-4, 8}]
```

`Power` also covers `scalar^Interval` and `Interval^Interval` (positive base):

```
In[4b]:= 2^Interval[{1, 3}]
Out[4b]= Interval[{2, 8}]

In[4c]:= Interval[{2, 3}]^Interval[{1, 2}]
Out[4c]= Interval[{2, 9}]
```

**Elementary functions.** `Exp`, `Log`, `Sqrt`, `Sinh`, `Tanh`, `ArcTan`,
`ArcSin`, `ArcSinh`, `ArcCosh`, `ArcTanh` (monotone increasing), `ArcCos`
(decreasing), and `Sin`, `Cos`, `Tan`, `Abs`, `Cosh` (non-monotone, with
critical-point range analysis) thread through intervals. The range analysis pins
an extremum to the exact value `±1` (or `Cosh`'s minimum `1`) when a
crest/trough lies inside the interval, and otherwise keeps the endpoint value
symbolic:

```
In[5]:= Sin[Interval[{2, 7}]]
Out[5]= Interval[{-1, Sin[2]}]

In[6]:= Sin[Interval[{2, 10}]]
Out[6]= Interval[{-1, 1}]

In[7]:= Sin[Interval[{2.5, 5.5}]]
Out[7]= Interval[{-1, 0.598472}]

In[8]:= Sqrt[Interval[{4, 9}]]
Out[8]= Interval[{2, 3}]
```

The threaded set also includes the **reciprocal trig/hyperbolic** functions
`Cot`, `Sec`, `Csc`, `Coth`, `Sech`, `Csch` (built as `1/base`, so a pole is
returned as a disjoint union), the **inverse-reciprocal** functions `ArcCot`,
`ArcSec`, `ArcCsc`, `ArcCoth`, `ArcSech`, `ArcCsch` (built as
`inverse_base(1/x)`, so a discontinuity such as `ArcCot`'s at 0 becomes a
disjoint union and an out-of-domain interval stays symbolic), the
**rounding/sign** functions `Sign`, `Floor`, `Ceiling` (monotone
non-decreasing), and several **special functions** on the sub-domain where each
is monotone: `Erf`/`Erfc` (all of ℝ), `Gamma`/`LogGamma`
(increasing above the minimum `x ≈ 1.4616`, decreasing below it), `Zeta` (on
`x > 1`), and `PolyGamma[n, ·]` (on `x > 0`, increasing for even `n`). Outside a
known monotone region a special function is left symbolic rather than risk a
wrong bound:

```
In[8b]:= Sec[Interval[{1, 2}]]
Out[8b]= Interval[{-Infinity, Sec[2]}, {Sec[1], Infinity}]

In[8c]:= Gamma[Interval[{2, 3}]]
Out[8c]= Interval[{1, 2}]

In[8d]:= Zeta[Interval[{2, 3}]]
Out[8d]= Interval[{Zeta[3], 1/6 Pi^2}]

In[8e]:= Gamma[Interval[{1, 2}]]        (* straddles the minimum: stays symbolic *)
Out[8e]= Gamma[Interval[{1, 2}]]
```

**General special/piecewise-function threading.** Beyond the hand-coded set
above, *any* `NumericFunction` applied to a single `Interval` argument threads by
**certified monotonicity**: the evaluator interval-evaluates the symbolic
derivative `D[f[x], x]` over each pair `[a, b]`, and if the enclosure is entirely
`≥ 0` (resp. `≤ 0`), `f` is monotone increasing (resp. decreasing) there and the
range is `[f(a), f(b)]` (swapped for decreasing). This is rigorous — interval
evaluation only ever over-estimates the derivative range, so a certified sign is
always correct, and an un-provable one falls back to symbolic. It reuses `D[]`
and the interval-aware evaluator, so no per-function monotonicity analysis is
needed:

```
In[8f]:= Erfi[Interval[{0, 1}]]
Out[8f]= Interval[{0, Erfi[1]}]

In[8g]:= ExpIntegralEi[Interval[{1, 2}]]
Out[8g]= Interval[{ExpIntegralEi[1], ExpIntegralEi[2]}]

In[8h]:= PolyLog[2, Interval[{0, 1/2}]]
Out[8h]= Interval[{0, -1/2 Log[2]^2 + 1/12 Pi^2}]
```

This covers `Erfi`, `ExpIntegralEi`, `LogIntegral`, `PolyLog[n, ·]`, the
incomplete `Gamma[a, ·]`, and any function whose derivative reduces to
already-threading heads. A few functions with a self-referential or awkward
derivative are handled instead by bespoke monotone-sub-domain rows: `InverseErf`
(on `(-1, 1)`), `InverseErfc` (on `(0, 2)`, decreasing), `ProductLog` (on
`(-1/e, ∞)`), and `HarmonicNumber` (on `(-1, ∞)`). The **piecewise/step**
functions `UnitStep`, `Ramp`, `Round`, `IntegerPart` (all non-decreasing, so
endpoint threading is a rigorous enclosure) join `Floor`/`Ceiling`.
Oscillatory heads (`BesselJ`/`BesselY`, `Sinc`, `SinIntegral`, `Fresnel*`, …) and
the sawtooth `FractionalPart` stay symbolic — no monotone bound is certifiable,
and the certifier terminates quickly rather than chasing a derivative chain that
never bottoms out.

**Endpoints and comparisons.** `Min[interval]` and `Max[interval]` return the
lowest and highest endpoints. Relational operators (`Equal`, `Less`, `Greater`,
…) yield explicit `True`/`False` when the interval is disjoint from the other
operand, and stay symbolic otherwise:

```
In[9]:= Min[Interval[{2, 7}]]
Out[9]= 2

In[10]:= Interval[{5, 8}] > Pi
Out[10]= True

In[11]:= Interval[{1, 4}] > Pi
Out[11]= Interval[{1, 4}] > Pi
```

Interval-vs-interval comparisons decide `True`/`False` when the two intervals
are disjoint, and stay symbolic when they overlap.

**Companion set operations.**

- `IntervalUnion[i1, i2, ...]` — the merged union (also `Interval`'s multi-argument
  form). `Flat` and `Orderless`.
- `IntervalIntersection[i1, i2, ...]` — the intersection, or the empty interval
  `Interval[]` when the arguments are disjoint. `Flat` and `Orderless`.
- `IntervalMemberQ[interval, x]` — `True` if the number `x` lies in `interval`;
  `IntervalMemberQ[interval, other]` tests whether the interval `other` is wholly
  contained.

```
In[12]:= IntervalUnion[Interval[{1, 3}], Interval[{2, 5}], Interval[{10, 11}]]
Out[12]= Interval[{1, 5}, {10, 11}]

In[13]:= IntervalIntersection[Interval[{1, 5}], Interval[{3, 8}]]
Out[13]= Interval[{3, 5}]

In[14]:= IntervalMemberQ[Interval[{1, 5}], Pi]
Out[14]= True
```

**Attributes.** `Interval` is `Protected`. It is deliberately *not*
`NumericFunction` (an interval is a set, not a scalar) and *not* `Listable`.

Implemented in `src/interval.{c,h}`; the arithmetic hooks live in
`src/plus.c` (`add_numbers`), `src/times.c` (the collector's interval
accumulator), and `src/power.c`.

**Not modelled.** Approximate numbers become *point* intervals
(`Interval[1.]` → `Interval[{1., 1.}]`) — a tight, rigorous enclosure of the
machine value; the reference's half-ULP "uncertainty preloading" of a machine number
is not replicated, because Mathilda's printer shows full precision and would
render `[1.−ulp, 1.+ulp]` as `{0.9999…, 1.0000…}` rather than `{1., 1.}`.
`Interval` is also not (yet) usable as a geometric `Region`, since Mathilda has
no `RegionMember`/region subsystem.
