# Elementary Functions

## Trig Functions
`Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc` and their inverses `ArcSin`, etc.

**Features**:
- `Listable`, `NumericFunction`.
- Exact values for rational multiples of `Pi` with denominators `1, 2, 3, 4, 5, 6, 10, 12`.
- `ArcTan[x, y]` computes the quadrant-aware inverse tangent.
- Numeric evaluation on `Complex[MPFR, MPFR]` (and out-of-real-domain
  MPFR arguments like `ArcSin[N[2, 50]]`) is carried at MPFR precision
  via the `numeric_mpfr_apply_complex_unary` helper rather than coerced
  to `double` through libc's `csin / ccos / ...` (see
  `src/numeric_complex.c`).

**Wrong arity.** Every elementary function validates its argument count and,
on a mismatch, emits the standard Mathematica diagnostic on stderr and leaves
the call unevaluated (returns the expression unchanged). The unary functions
(`Sin`, `Cos`, …, `ArcSin`, …, `Exp`, `Sqrt`, and the hyperbolic set below)
use `argx`; the variadic `ArcTan` and `Log` (1 or 2 arguments) use `argt`. The
message is produced by the shared `builtin_arg_error` helper in `src/common.c`.

```
Sin[]          -> Sin::argx: Sin called with 0 arguments; 1 argument is expected.
Sin[1, 2, 3]   -> Sin::argx: Sin called with 3 arguments; 1 argument is expected.
ArcTan[1,2,3]  -> ArcTan::argt: ArcTan called with 3 arguments; 1 or 2 arguments are expected.
```

```mathematica
In[1]:= Sin[Pi/6]
Out[1]= 1/2

In[2]:= ArcTan[1, 1]
Out[1]= 1/4*Pi
```

## Hyperbolic Functions
`Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`, `Csch` and their inverses.

**Features**:
- `Listable`, `NumericFunction`.
- Special values for `0` and `Infinity`.
- `Complex[MPFR, MPFR]` and out-of-real-domain MPFR inputs (e.g.
  `ArcCosh[N[0.5, 50]]`, `ArcTanh[N[2, 50]]`) preserve the working
  precision via the shared MPFR-complex helpers in
  `src/numeric_complex.c`.

## Exponential and Logarithmic Functions
- `Exp[z]`: Natural exponential. Canonicalizes `Exp[I*q*Pi]` (rational `q`) to `(-1)^q`, matching Mathematica — e.g. `Exp[I Pi/5] -> (-1)^(1/5)`, `Exp[I Pi] -> -1`, `Exp[I Pi/2] -> I`. `Power[-1, q]` reduces the half/integer cases itself and leaves irreducible roots intact instead of over-eagerly expanding into trig radicals.
- `Log[z]`: Natural logarithm. For a negative integer `n` (including `BigInt`), rewrites `Log[n]` as `I Pi + Log[-n]` (principal branch). For an exact pure-imaginary argument `Complex[0, b]` with `b` a real numeric, rewrites `Log[b I]` as `Log[Abs[b]] + Sign[b] (I Pi)/2` on the principal branch (e.g. `Log[I] = (I Pi)/2`, `Log[-3 I] = -((I Pi)/2) + Log[3]`); inexact imaginaries fall through to the numeric path. For a negative MPFR real or any `Complex[MPFR, MPFR]`, evaluates `log(hypot) + i atan2` at MPFR precision.
- `Log[b, z]`: Logarithm to base `b`.

**Zero arguments.** `Log` distinguishes exact zero (a directed limit) from
inexact zero (ambiguous direction in floating point), matching Mathematica:

- `Log[0]` -> `-Infinity`; `Log[b, 0]` -> `-Infinity` when `b` is a
  real-positive numeric with `b > 1` (or one of the known-positive
  symbols `E`, `Pi`), `+Infinity` when `0 < b < 1` (since `Log[b] < 0`).
- `Log[0.]` -> `Indeterminate`; `Log[b, 0.]` -> `Indeterminate` (both
  `+0.0` and `-0.0`).

**Wrong arity.** `Log` with 0, 3, or more arguments emits the standard
Mathematica diagnostic and leaves the call unevaluated:

```
Log[]            -> Log::argt: Log called with 0 arguments; 1 or 2 arguments are expected.
Log[2, 3, 4, 5]  -> Log::argt: Log called with 4 arguments; 1 or 2 arguments are expected.
```

```mathematica
In[1]:= Log[2, 8]
Out[1]= 3

In[2]:= Exp[I * Pi]
Out[2]= -1

In[3]:= Log[-5]
Out[3]= Log[5] + (I) Pi
```

**Power/Log cancellation.** `builtin_power` recognizes a `Log` in the
exponent and strips it:

- `Exp[c Log[a]]               -> a^c`   (since `Log[E] = 1`)
- `Power[base, c Log[base, a]] -> a^c`   (the `Log[base]` denominator cancels)

Internally `Log[b, a]` is represented as `Log[a] * Log[b]^(-1)`, so
`b^Log[b, a] -> a` and `Power[10, 3 Log[10, x]] -> x^3` both fall out of the
same rule. The coefficient `c` may be any (sub-)expression; it is the
product of whatever factors remain after removing the matching `Log[a]`
and (when `base != E`) `Log[base]^(-1)` factors from the exponent.

```mathematica
In[4]:= Exp[b Log[a]]
Out[4]= a^b

In[5]:= b^Log[b, a]
Out[5]= a
```

**Log of a Power with a matching base.** `builtin_log` folds
`Log[E^k] -> k` and `Log[b, b^k] -> k` when the exponent `k` is a real
numeric (integer, bigint, rational, or real) and (for the two-argument
form) `b` is a known-positive value (positive numeric, or a symbol
like `E` / `Pi`). The real-`k` guard keeps us on the principal branch --
for complex `k` the identity `Log[b^k] = k Log[b]` can fail by a
multiple of `2 Pi i`.

```mathematica
In[6]:= Log[E^4]
Out[6]= 4

In[7]:= Log[E^(1/3)]
Out[7]= 1/3

In[8]:= Log[2, 2^(1/3)]
Out[8]= 1/3
```

## Trig / inverse-trig and hyperbolic / inverse-hyperbolic cancellation

`builtin_sin`, `builtin_cos`, ..., `builtin_csc` and their hyperbolic
counterparts fold the direct forward-of-inverse identities:

- `Sin[ArcSin[x]]   -> x`,  `Cos[ArcCos[x]]   -> x`, ...
- `Sinh[ArcSinh[x]] -> x`,  `Cosh[ArcCosh[x]] -> x`, ...

These hold identically over the complex numbers because each `ArcF` is a
right inverse of `F` by construction. The opposite direction
(`ArcSin[Sin[x]]`, etc.) is *not* folded: it only reduces to `x` on each
function's principal domain.

The two-argument form `ArcTan[x, y]` is deliberately excluded from the
`Tan[ArcTan[...]]` rule -- `Tan[ArcTan[x, y]] = y/x`, not a single
variable -- so `Tan[ArcTan[3, 4]]` stays unevaluated.

```mathematica
In[9]:= Sin[ArcSin[x^2 + 1]]
Out[9]= 1 + x^2

In[10]:= Tanh[ArcTanh[z]]
Out[10]= z

In[11]:= ArcSin[Sin[x]]
Out[11]= ArcSin[Sin[x]]
```

### Singular special values

`ArcTanh` and `ArcCoth` have poles at `±1`, and Mathilda evaluates them to the
corresponding directed infinities (matching Mathematica):

```mathematica
In[12]:= ArcTanh[1]
Out[12]= Infinity

In[13]:= ArcTanh[-1]
Out[13]= -Infinity

In[14]:= ArcCoth[1]
Out[14]= Infinity
```

The negative cases are reached through the odd-function fold
(`ArcTanh[-1] -> -ArcTanh[1] -> -Infinity`). These closed-form values let the
definite integrator (`Integrate`, Newton–Leibniz) detect that an antiderivative
such as `-ArcTanh[x]` blows up at an interior pole, so a divergent integral like
`Integrate[1/(x^2-1), {x, -Infinity, Infinity}]` correctly warns
`Integrate::idiv` and returns unevaluated instead of a spurious finite value.

### Values at Infinity

The inverse-trig functions fold at `±Infinity` (Wolfram-faithful), so that both
top-level evaluation and `Limit` return closed forms:

```mathematica
In[15]:= ArcTan[Infinity]
Out[15]= Pi/2

In[16]:= ArcCot[-Infinity]
Out[16]= 0

In[17]:= ArcSec[Infinity]
Out[17]= Pi/2
```

`ArcTan[±Infinity] = ±Pi/2`, `ArcCot[±Infinity] = 0` (odd), `ArcSec[±Infinity] =
Pi/2`, `ArcCsc[Infinity] = 0`; negative arguments are reached through the
odd/`Pi`-minus fold. The circular functions `Sin`, `Cos`, `Tan`, … stay
unevaluated at real infinities (they oscillate). See `Limit` in the calculus
reference for how these feed the compose-at-infinity rule.

`ArcTan` also resolves *directionless* and *directed* infinities:
`ArcTan[ComplexInfinity] = Indeterminate` (the approach direction is unknown),
while `ArcTan[DirectedInfinity[d]] = ±Pi/2` picks the sign from the quadrant of
`d` — `Re[d] > 0 -> Pi/2`, `Re[d] < 0 -> -Pi/2`, and on the imaginary axis by
the sign of `Im[d]` (so `ArcTan[DirectedInfinity[I Sqrt[2]]] = Pi/2`). This is
what lets `Limit[2 ArcTan[Sqrt[(1 + x)/(1 - x)]], x -> 1]` reach `Pi` via the
compose-at-infinity rule.


## TrigToExp
Converts trigonometric and hyperbolic functions in `expr` to exponentials.
- `TrigToExp[expr]`

**Features**:
- `Listable`, `Protected`.
- Operates on both circular and hyperbolic functions, and their inverses.
- Rewrites recursively, including trig heads nested inside the *argument* of an
  outer trig head — e.g. `TrigToExp[Sin[Cos[x]]]` exponentializes the inner
  `Cos[x]` as well, matching Wolfram Language semantics.
- Automatically threads over lists, equations, inequalities, and logic functions.

```mathematica
In[1]:= TrigToExp[Cos[x]]
Out[1]= 1/2 E^(-I x) + 1/2 E^(I x)

In[2]:= TrigToExp[ArcSin[x]]
Out[2]= -I Log[I x + Sqrt[1 - x^2]]

In[3]:= TrigToExp[{Sin[x], Cos[x], Tan[x]}]
Out[3]= {(-1/2*I) E^((I) x) + (1/2*I) E^((-I) x), 1/2 E^((-I) x) + 1/2 E^((I) x), (-I) E^((I) x)/(E^((-I) x) + E^((I) x)) + (I) E^((-I) x)/(E^((-I) x) + E^((I) x))}
```

## ExpToTrig
Converts exponentials in `expr` to trigonometric or hyperbolic functions.
- `ExpToTrig[expr]`

**Features**:
- `Listable`, `Protected`.
- Tries when possible to give results that do not involve explicit complex numbers.
- ExpToTrig is natively the inverse of `TrigToExp`.
- Automatically threads over lists, equations, inequalities, and logic functions.

```mathematica
In[1]:= ExpToTrig[Exp[I x]]
Out[1]= Cos[x] + I Sin[x]

In[2]:= ExpToTrig[Log[1 + I x] - Log[1 - I x]]
Out[2]= 2 I ArcTan[x]

In[3]:= ExpToTrig[Exp[I x] == -1]
Out[3]= Cos[x] + I Sin[x] == -1
```

## TrigExpand
Expands trigonometric functions in `expr` by splitting up sums and integer
multiples that appear in arguments of circular and hyperbolic trigonometric
functions.
- `TrigExpand[expr]`

**Features**:
- `Listable`, `Protected`.
- Operates on both circular (`Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc`) and
  hyperbolic (`Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`, `Csch`) functions.
- Applies angle-addition formulas to `Sin[a + b + …]`, `Cos[a + b + …]`,
  `Sinh[a + b + …]`, `Cosh[a + b + …]` to a fixed point.
- Applies multiple-angle reductions to `Sin[n x]`, `Cos[n x]`, `Sinh[n x]`,
  `Cosh[n x]` for integer `n`, recursively reducing to `Sin[x]` / `Cos[x]` /
  `Sinh[x]` / `Cosh[x]`.
- `Tan`, `Cot`, `Sec`, `Csc` (and `Tanh`, `Coth`, `Sech`, `Csch`) with sum or
  integer-multiple arguments are rewritten as ratios of `Sin`/`Cos`
  (resp. `Sinh`/`Cosh`) and then expanded.
- Distributes products over sums via `Expand` so the result is a flat sum of
  monomials.
- Applies the Pythagorean identities `Sin[x]^2 + Cos[x]^2 -> 1` and
  `Cosh[x]^2 - Sinh[x]^2 -> 1` as a final reduction pass, including powers of
  both identities for any integer `n >= 1`:
    - `Sin[n x]^2 + Cos[n x]^2` expands to `(Sin[x]^2 + Cos[x]^2)^n` and
      collapses to `1` via a Factor-based reduction.
    - `Cosh[n x]^2 - Sinh[n x]^2` factors as
      `(Cosh[x] + Sinh[x])^n (Cosh[x] - Sinh[x])^n` and collapses to `1`.
  Negated and scalar-weighted forms (e.g. `-Sin[n x]^2 - Cos[n x]^2`,
  `-5 (Sin[n x]^2 + Cos[n x]^2)`, `Sinh[n x]^2 - Cosh[n x]^2`) collapse to the
  expected signed constant — the Pythagorean rules match both possible signs
  that `Factor` may emerge with and allow an arbitrary remainder of factors in
  the surrounding `Times`. Expressions that contain a denominator (any
  `Power[_, negative_Integer]` subterm) skip the Factor pass so that canonical
  forms such as `(2 Cos[x] Sin[x])/(Cos[x]^2 - Sin[x]^2)` are preserved.
  Inputs without a Pythagorean-eligible squared structure (no pair
  `Sin[a]^k`/`Cos[a]^k` or `Sinh[a]^k`/`Cosh[a]^k` with the same argument and
  `k >= 2`) likewise skip the Factor pass; the multivariate polynomials that
  multi-angle expansions such as `TrigExpand[Sin[2 x + 3 y]]` produce would
  otherwise make `Factor` prohibitively slow without yielding any collapse.
  The Factor pass is also skipped when the expanded form contains more than
  two distinct squared trigonometric atoms (e.g. `Cos[x]^2`, `Sin[x]^2`,
  `Cos[y]^2`, `Sin[y]^2` together): even if a Pythagorean pair is structurally
  present, `Factor` on the resulting dense multivariate polynomial stalls
  without producing a useful collapse.
- Automatically threads over lists (via `Listable`), as well as equations,
  inequalities (`Equal`, `Unequal`, `Less`, `LessEqual`, `Greater`,
  `GreaterEqual`, `SameQ`, `UnsameQ`), and logic functions (`And`, `Or`,
  `Not`, `Xor`, `Implies`).

```mathematica
In[1]:= TrigExpand[Sin[2 x]]
Out[1]= 2 Cos[x] Sin[x]

In[2]:= TrigExpand[Sin[x + y]]
Out[2]= Cos[x] Sin[y] + Cos[y] Sin[x]

In[3]:= TrigExpand[Sin[3 x]]
Out[3]= 3 Cos[x]^2 Sin[x] - Sin[x]^3

In[4]:= TrigExpand[Cos[x + y + z]]
Out[4]= -Cos[x] Sin[y] Sin[z] - Cos[y] Sin[x] Sin[z] - Cos[z] Sin[x] Sin[y] + Cos[x] Cos[y] Cos[z]

In[5]:= TrigExpand[Sin[x]^2 + Cos[x]^2]
Out[5]= 1

In[5b]:= TrigExpand[Sin[4 x]^2 + Cos[4 x]^2]
Out[5b]= 1

In[6]:= TrigExpand[Sinh[4 x]]
Out[6]= 4 Cosh[x] Sinh[x]^3 + 4 Cosh[x]^3 Sinh[x]

In[7]:= TrigExpand[Cosh[x - y]]
Out[7]= -Sinh[x] Sinh[y] + Cosh[x] Cosh[y]

In[8]:= TrigExpand[Tanh[2 t]]
Out[8]= (2 Cosh[t] Sinh[t])/(Cosh[t]^2 + Sinh[t]^2)

In[9]:= TrigExpand[{Tan[2 x], Sinh[x + y]}]
Out[9]= {(2 Cos[x] Sin[x])/(Cos[x]^2 - Sin[x]^2), Cosh[x] Sinh[y] + Cosh[y] Sinh[x]}

In[10]:= TrigExpand[1 < Cos[x + y] < 2]
Out[10]= 1 < -Sin[x] Sin[y] + Cos[x] Cos[y] < 2
```

## TrigFactor
Factors trigonometric functions in `expr`. Broadly the functional inverse of
`TrigExpand`: combines sums of trigonometric products into angle-sum forms,
collapses Pythagorean identities (both circular and hyperbolic), recognises
the reverse of the double-angle identities, and factors polynomial structure
over the trigonometric atoms.
- `TrigFactor[expr]`

**Features**:
- `Listable`, `Protected`.
- Operates on both circular (`Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc`) and
  hyperbolic (`Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`, `Csch`) functions.
- Pipeline:
  1. Rewrite reciprocal heads (`Tan`, `Cot`, `Sec`, `Csc`, and their
     hyperbolic analogs) as `Sin`/`Cos`/`Sinh`/`Cosh` ratios so that `Factor`
     sees the full polynomial structure.
  2. Combine into a single rational via `Together`.
  3. Run `Factor` on the resulting rational; trigonometric atoms are treated
     as independent polynomial variables. The `Factor` pass is skipped when
     the post-`Together` form contains more than two distinct squared
     trigonometric atoms (e.g. `Sin[x]^2`, `Cos[x]^2`, `Sinh[y]^2`,
     `Cosh[y]^2` together): on such dense multivariate polynomials Factor's
     trial-division loop stalls without producing a useful factorization, and
     the identity rules in step 4 still match Pythagorean structure that
     survives in the post-`Together` factored form (e.g.
     `(Sin[x]^2 + Cos[x]^2)(Cosh[y]^2 - Sinh[y]^2)` collapses directly to
     `1` via the `Times`-context Pythagorean rules).
  4. Apply identity collapse rules via `ReplaceRepeated`: Pythagorean
     identities (`Sin^2 + Cos^2 -> 1`, `Cosh^2 - Sinh^2 -> 1`, with and
     without arbitrary coefficients), reverse angle-addition
     (`Sin[a]Cos[b] ± Cos[a]Sin[b] -> Sin[a ± b]`,
     `Cos[a]Cos[b] ± Sin[a]Sin[b] -> Cos[a ∓ b]`, and hyperbolic analogs),
     reverse double-angle (`2 Sin Cos -> Sin[2x]`,
     `Cos^2 - Sin^2 -> Cos[2x]`, `Cosh^2 + Sinh^2 -> Cosh[2x]`), and
     factored-form variants such as `(Cos - Sin)(Cos + Sin) -> Cos[2x]`,
     `(Cosh - 1)(Cosh + 1) -> Sinh^2`, and
     `(Cosh - Sinh)(Cosh + Sinh) -> 1` that arise naturally from `Factor`.
  5. Restore `Tan`/`Cot`/`Sec`/`Csc` (and hyperbolic analogs) from the
     `Sin`/`Cos` ratio form so reciprocal heads survive the round-trip.
- Two paths are tried: the primary pipeline (preserves angle-sum structure)
  and a fallback that `TrigExpand`s the argument first (catches
  cancellations that only become visible after the angle-sum is expanded,
  e.g. `Cos[x + y] + Sin[x] Sin[y] -> Cos[x] Cos[y]`). The fallback runs
  only when the primary pipeline leaves the expression unchanged, so
  structurally productive inputs (e.g. `Sin[x + y]^2 + Tan[x + y]`) avoid
  the expensive expanded-rational path. The final result is the smaller of
  the two by leaf count; ties favour the primary pipeline.
- Automatically threads over lists (via `Listable`), as well as equations,
  inequalities (`Equal`, `Unequal`, `Less`, `LessEqual`, `Greater`,
  `GreaterEqual`, `SameQ`, `UnsameQ`), and logic functions (`And`, `Or`,
  `Not`, `Xor`, `Implies`).

```mathematica
In[1]:= TrigFactor[Sin[x]^2 + Cos[x]^2]
Out[1]= 1

In[2]:= TrigFactor[Cosh[x]^2 - Sinh[x]^2]
Out[2]= 1

In[3]:= TrigFactor[2 Sin[x] Cos[x]]
Out[3]= Sin[2 x]

In[4]:= TrigFactor[Cos[x]^2 - Sin[x]^2]
Out[4]= Cos[2 x]

In[5]:= TrigFactor[Sin[a] Cos[b] + Cos[a] Sin[b]]
Out[5]= Sin[a + b]

In[6]:= TrigFactor[Cos[a] Cos[b] + Sin[a] Sin[b]]
Out[6]= Cos[a - b]

In[7]:= TrigFactor[Sin[x]^2 + Tan[x]^2]
Out[7]= Tan[x]^2 (1 + Cos[x]^2)

In[8]:= TrigFactor[Cosh[x]^2 - Cosh[x]^4]
Out[8]= -Cosh[x]^2 Sinh[x]^2

In[9]:= TrigFactor[Sin[x+y]^2 + Tan[x+y]]
Out[9]= Tan[x + y] (1 + Cos[x + y] Sin[x + y])

In[10]:= TrigFactor[Cos[x + y] + Sin[x] Sin[y]]
Out[10]= Cos[x] Cos[y]

In[11]:= TrigFactor[Cos[x]^4 - Sin[x]^4]
Out[11]= Cos[2 x]

In[12]:= TrigFactor[{Sin[x]^2 + Cos[x]^2, 2 Sinh[x] Cosh[x]}]
Out[12]= {1, Sinh[2 x]}

In[13]:= TrigFactor[Sin[x]^2 + Cos[x]^2 == 1]
Out[13]= True
```

## TrigReduce

`TrigReduce[expr]` rewrites products and integer powers of single-argument
trigonometric functions in `expr` in terms of trigonometric functions with
combined arguments. It is the inverse-of-`TrigExpand` direction for the
product/power identities: given a trigonometric polynomial, `TrigReduce`
typically yields a linear expression involving trigonometric functions
with more complicated arguments.

- `TrigReduce[expr]`

**Features**:
- Applies the classical product-to-sum identities (Sin·Cos, Sin·Sin,
  Cos·Cos, plus the four hyperbolic analogues) and the power-reduction
  identities (Sin² → (1 − Cos[2x])/2, etc., extended recursively to
  any positive integer power).
- Operates on both circular and hyperbolic functions; rewrites
  `Tan`/`Cot`/`Sec`/`Csc` (and the hyperbolic reciprocals) as `Sin`/`Cos`
  ratios before reduction, then restores the reciprocal head where the
  result has the matching shape.
- Recognises angle-addition forms produced after `Together`
  (`Sin[a]·Cos[b] + Cos[a]·Sin[b] → Sin[a + b]` and analogues), with
  sign variants and hyperbolic counterparts.
- Includes a sign-cancellation pass for the `Sin[a − b] + Sin[b − a] = 0`
  shape that arises when the product-to-sum rules bind asymmetrically;
  the same pass handles the corresponding `Cos`/`Sinh`/`Cosh` parities.
- Applies an `Expand`/`Together` cleanup at the end so trivial outer
  fractions (`1/2 (2 X + 2 Y)`) flatten while genuine rationals
  (`(3 − 4 Cos[2 x] + Cos[4 x])/2`) survive in normalised form.
- Memoised through the same `FactorMemo` mechanism used by
  `TrigExpand`/`TrigFactor`/`TrigToExp`, so repeated invocations during
  `Simplify` candidate-set search amortise.
- `Listable`, plus explicit threading over `Equal`, `Unequal`, `Less`,
  `LessEqual`, `Greater`, `GreaterEqual`, `SameQ`, `UnsameQ`, `And`,
  `Or`, `Not`, `Xor`, `Implies` (mirrors `TrigExpand` / `TrigFactor`).
- Idempotent on already-reduced inputs and a no-op on non-trig
  expressions or single trig calls of compound arguments.

```mathematica
In[1]:= TrigReduce[2 Cos[x]^2]
Out[1]= 1 + Cos[2 x]

In[2]:= TrigReduce[2 Sin[x] Cos[y]]
Out[2]= Sin[x + y] + Sin[x - y]

In[3]:= TrigReduce[2 Cosh[x] Cosh[y]]
Out[3]= Cosh[x + y] + Cosh[x - y]

In[4]:= TrigReduce[Sin[a] (Cos[b] - Sin[b]) + Cos[a] (Sin[b] + Cos[b])]
Out[4]= Cos[a + b] + Sin[a + b]

In[5]:= TrigReduce[Tan[x] + Tan[y]]
Out[5]= Sec[x] Sec[y] Sin[x + y]

In[6]:= TrigReduce[Coth[x] + Coth[y]]
Out[6]= Csch[x] Csch[y] Sinh[x + y]

In[7]:= TrigReduce[Sin[x]^4]
Out[7]= 1/8 (3 - 4 Cos[2 x] + Cos[4 x])

In[8]:= TrigReduce[2 Sin[x + y] Cos[x - y]]
Out[8]= Sin[2 x] + Sin[2 y]

In[9]:= TrigReduce[{Tan[x] + Cot[y], Tanh[x] - Coth[y]}]
Out[9]= {Cos[x - y] Csc[y] Sec[x], -Cosh[x - y] Csch[y] Sech[x]}

In[10]:= TrigReduce[4 Sin[x]^4 == 1 && 2 Cos[x]^2 >= 1]
Out[10]= 1/2 (3 - 4 Cos[2 x] + Cos[4 x]) == 1 && 1 + Cos[2 x] >= 1
```

`TrigReduce` is also offered as a `Simplify` transform: when the search
sees an angle-addition expansion such as
`Sin[a] (Cos[b] − Sin[b]) + Cos[a] (Sin[b] + Cos[b])`, the reduced form
`Cos[a + b] + Sin[a + b]` has fewer leaves and the score-gate selects
it.

## Piecewise and Rounding Functions
`Floor`, `Ceiling`, `Round`, `IntegerPart`, `FractionalPart`.

**Features**:
- `Listable`.
- Applied component-wise to `Complex` numbers.
- `Round` rounds to the nearest even integer for ties.
- `Floor[x, a]` returns the greatest multiple of `a` $\le x$.
- **Exact real numeric arguments** that no leaf branch resolves (e.g. `Pi`, `E`, surds, and products such as `10000000 3^(2/3)` or `25000000000000000000 Pi`) are numericalized to MPFR at increasing precision and reduced to the exact integer. The result is only accepted once two successive precisions (starting at 256 bits, doubling to a 65536-bit cap) agree on the integer — an interval-style certification that never returns a wrong answer; a value that cannot be certified (or is not a pure real number) is left unevaluated. For `FractionalPart` the answer is kept exact as `x - IntegerPart[x]` (matching Mathematica).

**Symbolic simplifications** (`Floor`, `Ceiling`, `Round` only -- `IntegerPart` and `FractionalPart` are excluded):

*Sign extraction* -- pulls a leading `-1` out through the rounding head, swapping `Floor`/`Ceiling` (and leaving `Round` unchanged):
- `Floor[-x]   :> -Ceiling[x]`
- `Ceiling[-x] :> -Floor[x]`
- `Round[-x]   :> -Round[x]`

The trigger is `expr_is_superficially_negative`, so `Floor[-2 x]` and `Ceiling[-3 x y]` reduce as well.

*Idempotency / composition* -- the inner expression is already an integer, so the outer rounding is a no-op:
- `Floor[Floor[x]]   :> Floor[x]`
- `Ceiling[Ceiling[x]] :> Ceiling[x]`
- `Round[Round[x]]   :> Round[x]`
- `Floor[Ceiling[x]] :> Ceiling[x]`
- `Ceiling[Floor[x]] :> Floor[x]`
- `Floor[Round[x]]   :> Round[x]`
- `Ceiling[Round[x]] :> Round[x]`
- `Round[Floor[x]]   :> Floor[x]`
- `Round[Ceiling[x]] :> Ceiling[x]`

These rules compose under fixed-point evaluation, so e.g. `Ceiling[Floor[Ceiling[x]]] -> Ceiling[x]`, and `Floor[-x] + Ceiling[x] -> 0`.

```mathematica
In[1]:= Round[2.5]
Out[1]= 2

In[2]:= Round[3.5]
Out[2]= 4

In[3]:= Floor[-x] + Ceiling[x]
Out[3]= 0

In[4]:= Ceiling[Floor[Ceiling[x]]]
Out[4]= Ceiling[x]

In[5]:= Round[10000000*3^(2/3)]
Out[5]= 20800838

In[6]:= Round[25000000000000000000 Pi]
Out[6]= 78539816339744830962

In[7]:= FractionalPart[10000000*3^(2/3)]
Out[7]= -20800838 + 10000000 3^(2/3)
```

## UnitStep

`UnitStep[x]` is the unit step (Heaviside) function: `0` for $x < 0$ and `1`
for $x \ge 0$. The value at zero is `1`.

`UnitStep[x1, x2, ...]` is the multidimensional unit step, equal to `1` only
when none of the `xi` are negative, and `0` otherwise (it behaves as the
product of the per-argument `UnitStep[xi]`). `UnitStep[]` is `1`.

**Features**:
- `Listable`, `NumericFunction`, `Orderless`, `Protected`.
- The result is **always exact** -- an integer `0` or `1` -- for real numeric
  input, including `Real`/`MPFR` arguments (e.g. `UnitStep[{-1.6, 3.2}]` gives
  `{0, 1}`).
- **Exact symbolic real arguments** (`Pi`, `Sqrt[2]`, `E - 3`, ...) are
  resolved by numerical certification: the argument is numericalized to MPFR at
  increasing precision and the sign is accepted only once two successive
  precisions agree on the same non-zero sign. This separates tight cases such
  as `Sqrt[2] - 99/70` ($\approx -6.4\times10^{-5}$) from zero without ever
  guessing; an argument whose sign cannot be certified is left unevaluated.
- Non-real arguments (a `Complex` with non-zero imaginary part) and unresolved
  symbolic arguments are left unevaluated. In a multidimensional call the
  proven-non-negative arguments are dropped (they contribute a factor of `1`),
  so e.g. `UnitStep[1, x]` reduces to `UnitStep[x]`.

**Derivative** -- via the product rule, each argument contributes
`Piecewise[{{Indeterminate, xi == 0}}, 0]`:

```mathematica
In[1]:= UnitStep[0]
Out[1]= 1

In[2]:= UnitStep[1, Pi, 5.3]
Out[2]= 1

In[3]:= UnitStep[{-1.6, 3.200000000000}]
Out[3]= {0, 1}

In[4]:= UnitStep[Sqrt[2] - 99/70]
Out[4]= 0

In[5]:= D[UnitStep[x], x]
Out[5]= Piecewise[{{Indeterminate, x == 0}}, 0]

In[6]:= D[UnitStep[x, y, z], z]
Out[6]= UnitStep[x, y] Piecewise[{{Indeterminate, z == 0}}, 0]
```

## Chop

`Chop[expr]` replaces approximate real numbers in `expr` whose absolute
value is below `10^-10` by the exact integer `0`.

`Chop[expr, delta]` uses `|delta|` as the threshold instead of the
default.

**Features**:
- `Protected`.
- Walks the entire expression tree, so small real-valued subterms inside
  arbitrary heads, lists, and held forms are all chopped.
- Exact numbers (`Integer`, BigInt, `Rational`, symbolic constants) and
  symbolic input pass through untouched.
- `delta` may be supplied as `Integer`, `Real`, `Rational[n, d]`,
  BigInt, or (when `USE_MPFR=1`) `MPFR`; its absolute value is the
  effective tolerance.

**Complex handling**.  `Complex[re, im]` whose components are both
machine reals is the "machine complex" case and gets Mathematica's
special treatment:

| `re` below tolerance | `im` below tolerance | result |
|----------------------|----------------------|--------|
| yes                  | yes                  | exact integer `0` |
| no                   | yes                  | `re` (Complex wrapper dropped, machine real survives) |
| yes                  | no                   | `Complex[0., im]` &mdash; the real part is the machine zero `0.0`, not an exact `0`, preserving the machine-complex shape |
| no                   | no                   | unchanged |

When `Complex[re, im]` has at least one exact component (e.g.
`Complex[1, 1.*^-12]`), Chop recurses into each part normally: a tiny
`Real` becomes the exact integer `0`, and `builtin_complex` then
collapses `Complex[r, 0]` to `r` on the next evaluator pass.

```mathematica
In[1]:= Chop[Exp[N[Range[4] Pi I]]]
Out[1]= {-1., 1., -1., 1.}

In[2]:= Chop[N[Pi] - Rationalize[N[Pi], 10^-12]] === 0
Out[2]= True

In[3]:= Chop[N[Pi] - Rationalize[N[Pi], 10^-12], 10^-14] === 0
Out[3]= False

In[4]:= Chop[10.^-12 + 2. I]
Out[4]= 0. + 2. I

In[5]:= Chop[2. + 10.^-12 I]
Out[5]= 2.
```

## Clip

`Clip[x]` clamps a numeric value to the closed interval `[-1, 1]`.

`Clip[x, {min, max}]` clamps to `[min, max]`.

`Clip[x, {min, max}, {vmin, vmax}]` returns `vmin` when `x < min`,
`vmax` when `x > max`, and `x` otherwise.  The replacement values
need not be numeric.

**Features**:
- `NumericFunction`, `Protected`.
- Threads over a `List` in the first argument: `Clip[{x1, x2, ...}, ...]`
  maps Clip element-wise over the list.  The `{min, max}` and
  `{vmin, vmax}` configuration lists are explicitly **not** threaded
  over -- threading is implemented inside the builtin (not via the
  `Listable` attribute) so the bounds and replacement lists stay intact.
- Symbolic numeric constants (`Pi`, `E`, etc.) are numericalized via
  `N` only to decide which side of the interval `x` lies on; the
  original symbolic `x` is returned when `min <= x <= max`, never
  the numeric approximation.
- `Infinity` and `-Infinity` are handled directly: `Clip[Infinity]`
  yields `vmax` (or the default `1`), `Clip[-Infinity]` yields `vmin`.
- Complex (non-real) input emits a one-shot `Clip::ncompl` warning and
  the call stays unevaluated.  Use `Re[z]`, `Im[z]` to clip the parts
  separately.
- Symbolic input for which the position cannot be decided
  numerically (e.g. `Clip[a]`) stays unevaluated so user-supplied
  rules can intercept it.

```mathematica
In[1]:= Clip[7.5]
Out[1]= 1

In[2]:= Clip[-5/2, {-2, 2}]
Out[2]= -2

In[3]:= Clip[Pi, {-9, 7}, {11, 28}]
Out[3]= Pi

In[4]:= Clip[{-2, 0, 2}]
Out[4]= {-1, 0, 1}

In[5]:= Clip[Infinity]
Out[5]= 1

In[6]:= Clip[2 - 3 I]
Clip::ncompl: Symbolic or noncomplex numerical arguments are expected.
Out[6]= Clip[2 - 3 I]

In[7]:= Clip[Re[2 - 3 I]] + Clip[Im[2 - 3 I]] I
Out[7]= 1 - I

In[8]:= N[Clip[1/11, {1/7, 5}], 50]
Out[8]= 0.14285714285714285714285714285714285714285714285714
```

