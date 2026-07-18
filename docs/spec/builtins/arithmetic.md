# Arithmetic

Core arithmetic: the basic operators (`Plus`, `Times`, `Power`, `Sqrt`), exact and approximate numbers (`N`, precision control), combinatorial functions (`Factorial`, `Binomial`, `Fibonacci`), digit and radix manipulation, and complex-number component functions. Number-theoretic functions (`GCD`, `PowerMod`, `PrimeQ`, `FactorInteger`, ...) live in [`number-theory.md`](number-theory.md).

## Numeric literals (parser-level)


Numeric literals are classified by the parser based on the precision
the user typed.

- **Integer / BigInt** — literals with no `.`, `e`, `E`, or `*^`.
  Digit-count drives the choice between `EXPR_INTEGER` (`int64_t`) and
  `EXPR_BIGINT` (GMP).
- **`EXPR_REAL`** — real literals whose **implied precision** is at
  most `MachinePrecision` (≈ 15.95 decimal digits).
- **`EXPR_MPFR`** (arbitrary precision) — real literals whose implied
  precision exceeds `MachinePrecision`. The MPFR bit count is
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

Examples (verified against Mathematica):

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

Mathematica-style backtick suffixes (`` ` ``, `` `` ``) still override
the implicit decision: `3.14`50` is 50-digit MPFR, `3.14``49` is
49-digit accuracy, and `3``50` is `EXPR_REAL`.

## N

Numerical evaluation.
- `N[expr]` — machine-precision approximation.
- `N[expr, p]` — `p`-digit (arbitrary-precision, MPFR) approximation.

**Precision semantics**:
- `N[expr]` numericalizes only the **exact** parts of `expr` (integers,
  rationals, constants such as `Pi`, `E`, `GoldenRatio`, …) and **leaves the
  precision of numbers that are already approximate untouched**. So
  `N[N[Pi, 100]]` (i.e. `N[Pi, 100] // N`) stays at 100 digits rather than
  collapsing to `MachinePrecision`, matching Mathematica. Likewise `N[2.5`100]`
  keeps its 100-digit precision.
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
- Rationals whose numerator or denominator overflow `int64`/`double` still
  numericalize to a **single** real, not a frozen `Rational[Real, Real]`
  (`N[1/10^30]` → `1.e-30`, `N[10^400/3]` → `3.33e+399`). Out-of-range machine
  results are carried as a machine-precision MPFR real so the arbitrary
  exponent survives (matching `N[1001!]`).

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
- Returns `0` if no arguments are provided.
- Returns `Overflow[]` if integer addition overflows or if any argument is `Overflow[]`.

```mathematica
In[1]:= 1 + 2 + x + 2*x
Out[1]= 3 + 3*x
```

## Total

Gives the total of elements in a list.
- `Total[list]`
- `Total[list, n]`
- `Total[list, {n}]`
- `Total[list, {n1, n2}]`

**Features**:
- `Protected`.
- `Total[list]` is equivalent to `Apply[Plus, list]`.
- `Total[list, n]` totals all elements down to level `n`.
- `Total[list, {n}]` totals elements at level `n` only.
- Supports negative levels to count from the bottom (`-1` is the last dimension).
- Handles ragged arrays correctly by summing from the inside out when multiple levels are specified.
- `Total[list, Infinity]` totals all atoms in the expression.

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

## Times (*)

Symbolic product.
- `a * b * ...`

**Features**:
- `Flat`, `Orderless`, `Listable`.
- Combines numeric constants and groups identical bases into `Power` expressions.
- Handles `I` as `Complex[0, 1]`.
- Returns `1` if no arguments are provided.
- Returns `Overflow[]` if integer multiplication overflows or if any argument is `Overflow[]`.
- **Sqrt-coefficient absorption**: a non-trivial rational/integer coefficient
  combined with a single `Power[r, -1/2]` group folds into the canonical
  Mathematica form `sign(c) * Sqrt[c^2 / r]`, then Power's existing
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

## Power (^)

Exponentiation.
- `base ^ exp`

**Features**:
- `Listable`.
- Simplifies integer powers of integers.
- Returns `Overflow[]` if the result exceeds 64-bit integer limits.
- Reduces radicals (e.g., `8^(1/2)` becomes `2*Sqrt[2]`).
- Supports complex results for negative bases (e.g., `(-1)^(1/2)` becomes `I`).
  Higher-power cases for `q == 2` now also reduce: `(-1)^(3/2) → -I`,
  `(-1)^(5/2) → I`, `(-12)^(3/2) → -24 I Sqrt[3]` (the principal-branch
  rule `(-n)^(p/2) = I^p · |n|^(p/2)`).
- For a negative base, the residual `(-1)^(b/q)` exponent is reduced into
  `[0, 1)` by **floor** division of `p/q`, pulling out a `(-1)^a = ±1` sign
  that merges into the coefficient (Mathematica canonical form). This covers
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
- For `Power[Integer, Rational]` with positive base and positive
  `p/q` exponent, splits the base's prime factorisation into a
  product of distinct-prime powers grouped by reduced effective
  exponent (Mathematica canonical form). Triggers only when the
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

## Mod, Quotient, QuotientRemainder

- `Mod[n, m]`: Remainder of `n/m`.
- `Quotient[n, m]`: Integer part of `n/m`.
- `QuotientRemainder[n, m]`: Returns `{Quotient[n, m], Mod[n, m]}`.

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
- Works seamlessly for both machine integers and arbitrary-precision
  bignums (digits are computed in GMP and demoted to machine ints
  whenever they fit).

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
diagnostic and echo the call back unevaluated, matching Mathematica:

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
diagnostic and echo the call back unevaluated, matching Mathematica:

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
  digit `0` is **last**, matching Mathematica's convention.
- `DigitCount[n, b, d]`: Scalar count of digit `d` (an integer with
  `0 <= d < b`) in the base-`b` representation of `n`.

**Features**:
- `Protected` (intentionally not `Listable`). `DigitCount[{1,2,3}]` is left
  unevaluated rather than threading, because the natural result of
  `DigitCount[n]` is itself a list.
- `DigitCount[0]` is a list of zeros of length `b` (zero has no significant
  digits). In the 3-argument form, `DigitCount[0, b, d] == 0`.
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

**Features**:
- `Protected` (intentionally not `Listable`: the first argument *is* a list).
- Inverse of `IntegerDigits` / `IntegerString`. Since `IntegerDigits` discards
  the sign of `n`, `FromDigits[IntegerDigits[n]]` is `Abs[n]`, not `n`.
- Digits in the list (and characters in the string) need *not* be smaller
  than the base; they are carried via Horner's evaluation, matching
  Mathematica (e.g. `FromDigits[{7, 11, 0, 0, 0, 122}] == 810122`,
  `FromDigits["1A3C"] == 2042`).
- The all-integer (digits + base) case is computed exactly in GMP and
  demoted to a machine integer when it fits in `int64`; arbitrarily large
  bignums are returned otherwise.
- Symbolic, Real, or Rational digits or base trigger the polynomial Horner
  expansion `d[0] b^(n-1) + d[1] b^(n-2) + ... + d[n-1]`, which the
  evaluator simplifies normally. This single code path handles symbolic
  bases, symbolic digits, and inexact bases uniformly. For a string input,
  the base is required to be an integer; symbolic bases over strings leave
  the call unevaluated (silently).
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
  `len` least-significant digits are returned (matching Mathematica).

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
diagnostic and echo the call back unevaluated, matching Mathematica:

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
  the input precisions.  This matches Mathematica's contagion: an
  explicit `N[..., p]` lifts the exponent to the same precision.
- The base must be a real number `> 1`; non-positive, `<= 1`, or complex
  bases emit `RealExponent::ibase` and leave the call unevaluated.
- Complex arguments with non-zero imaginary part emit
  `RealExponent::realx` and leave the call unevaluated.
- Sign of `x` is discarded.
- Zero handling (Mathematica-compatible):
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

## Binomial

Gives the binomial coefficient $\binom{n}{m}$, generalised via
$\Gamma(n+1)/(\Gamma(m+1)\,\Gamma(n-m+1))$.

- `Binomial[n, m]`

**Features**:
- `Protected`, `Listable`, `NumericFunction`.
- Exact integer/integer path uses GMP (`mpz_bin_ui`), including the
  Pascal extension `Binomial[n, m] = (-1)^m Binomial[m-n-1, m]` for
  negative `n`.
- Machine-precision real branch via `tgamma`.
- Symmetric identity: when `n - m` simplifies to a non-negative
  integer `k ≤ 32`, reduces to `Binomial[n, k]` and expands as a
  falling-factorial polynomial. This catches `Binomial[n, n - 1] → n`,
  `Binomial[9/2, 7/2] → 9/2`, `Binomial[n + 1, n - 1] → n (n + 1)/2`,
  etc.
- Concrete non-negative integer `m ≤ 32` with any other `n` (symbolic,
  rational, complex, …) expands to the falling-factorial polynomial
  `n (n-1) (n-2) ... (n-m+1) / m!`, which the `Times`/`Plus` folders
  then simplify.

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
  precision.

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
