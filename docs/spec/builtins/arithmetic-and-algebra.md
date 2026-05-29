# Arithmetic and Algebra

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
- `Power[-1, p/q]` with even `q ≥ 4` and `|p| ≥ q` does integer-part
  extraction, e.g. `(-1)^(5/4) → -(-1)^(1/4)`, `(-1)^(7/6) → -(-1)^(1/6)`.
  `(-1)^(p/q)` with `|p| < q` (and odd-`q` negative-base cases for any
  base) continue to canonicalise via the existing path.  Negative bases
  other than `-1` with even `q ≥ 4` (e.g. `(-16)^(1/4)`) are still left
  unevaluated.
- Distributes power over product if the exponent is an integer.
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

## Numerator
Gives the numerator of an expression.
- `Numerator[expr]`

**Features**:
- `Protected`, `Listable`.
- Picks out terms which do not have superficially negative exponents.
- Can be used on rational and complex numbers.

```mathematica
In[1]:= Numerator[(x-1)(x-2)/(x-3)^2]
Out[1]= (-2 + x) (-1 + x)

In[2]:= Numerator[3/7 + I/11]
Out[2]= 33 + 7 I
```

## Denominator
Gives the denominator of an expression.
- `Denominator[expr]`

**Features**:
- `Protected`, `Listable`.
- Picks out terms which have superficially negative exponents.
- Can be used on rational and complex numbers.

```mathematica
In[1]:= Denominator[(x-1)(x-2)/(x-3)^2]
Out[1]= (-3 + x)^2

In[2]:= Denominator[3/7 + I/11]
Out[2]= 77
```

## Cancel
Cancels out common factors in the numerator and denominator of an expression.
- `Cancel[expr]`
- `Cancel[expr, Extension -> alpha]`

**Features**:
- `Protected`, `Listable`.
- Threads over equations, inequalities, logic functions, and sums dynamically.
- Evaluates greatest common divisors via polynomial GCD derivations avoiding extraneous expansions.
- Handles a single symbolic base appearing with rational fractional exponents (e.g. `Sqrt[y]`, `y^(1/3)`) by treating it as an algebraic generator: substitutes `y -> g^m` where `m` is the LCM of denominators, runs the polynomial cancellation in `g`, then substitutes back.
- The algebraic-generator pass runs `Together` on the substituted form (not just GCD-cancellation), so inputs whose `g`-substituted denominator is a Plus of terms with different `g`-denominators (e.g. `1/(g^2 - 1/g)` from `1/(y^(2/3) - 1/y^(1/3))`) are handled correctly.
- Extracts algebraic-constant atoms (`Sqrt[2]`, `Sqrt[3]`, `CubeRoot[5]`, `2^(2/3)`, ...) that appear in every summand of *both* numerator *and* denominator and divides them out before the polynomial GCD step. Closes a long-standing gap where `PolynomialGCD[Sqrt[2], Sqrt[2] + Sqrt[2] x^4]` returns 1 (the integer-content recursion treats `Sqrt[2]` as having content 1), so cancellations whose only shared factor was an algebraic constant survived as-is. The pass is intentionally narrow — only `Power[integer, rational/non-integer]` factors are eligible — to avoid disturbing the rational-function intermediates that the integration dispatcher pattern-matches against.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan) cancels common factors over `Q(alpha)` instead of `Q`. Implementation: lifts numerator and denominator into `Q(alpha)[x]` via the QAUPoly machinery, runs `qaupoly_gcd`, divides both sides by `g`, and re-renders. Works for single-fraction inputs; `Plus` inputs (sums of fractions) currently fall back to the no-extension path because `PolynomialQuotient` does not yet accept `Extension` (Phase 0.5 follow-up).

```mathematica
In[1]:= Cancel[(x^2 - 1) / (x - 1)]
Out[1]= 1 + x

In[2]:= Cancel[(x - y)/(x^2 - y^2) + (x^3 - 27)/(x^2 - 9)]
Out[2]= (9 + 3 x + x^2)/(3 + x) + 1/(x + y)

In[3]:= Cancel[(y - 1)/(Sqrt[y] - 1)]
Out[3]= 1 + Sqrt[y]

In[4]:= Cancel[(y - 1)/(y^(1/3) - 1)]
Out[4]= 1 + y^(1/3) + y^(2/3)

In[5]:= Cancel[1/(y^(2/3) - 1/y^(1/3))]
Out[5]= y^(1/3)/(-1 + y)

In[6]:= Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Sqrt[2]]
Out[6]= Sqrt[2] + x

In[7]:= Cancel[(x^3 - 2)/(x - 2^(1/3)), Extension -> 2^(1/3)]
Out[7]= 2^(2/3) + 2^(1/3) x + x^2

In[8]:= Cancel[Sqrt[2]/(Sqrt[2] + Sqrt[2] x^4)]
Out[8]= 1/(1 + x^4)

In[9]:= Cancel[(Sqrt[2] x + Sqrt[2] y)/(Sqrt[2] + Sqrt[2] x)]
Out[9]= (x + y)/(1 + x)
```

## Together
Puts terms in a sum over a common denominator, and cancels factors in the result.
- `Together[expr]`
- `Together[expr, Extension -> alpha]`

**Features**:
- `Protected`, `Listable`.
- Makes a sum of terms into a single rational function.
- Computes lowest common multiples (LCM) of denominators securely without unconditionally destroying pre-factored bases unnecessarily.
- Handles a single symbolic base appearing with rational fractional exponents (e.g. `y^(1/3)`, `y^(2/3)`, `y^(73/24)`) by treating it as an algebraic generator: substitutes `y -> g^m` where `m` is the LCM of denominators, runs the polynomial pipeline in `g`, then substitutes back.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan) combines into a single fraction with the standard combiner, then runs `Cancel[..., Extension -> alpha]` on the result so algebraic-coefficient cancellations fire. Effective on simple inputs like `1/(x - Sqrt[2]) + 1/(x + Sqrt[2])`, which collapses to `(2 x)/(x^2 - 2)`. Inputs whose summands themselves carry algebraic-coefficient denominators are deferred to Phase 0.5 (which will plumb `Extension` through `PolynomialLCM` / `PolynomialQuotient` / `together_recursive`).
- **Option `Extension -> Automatic` with polynomial radicands** (Phase E, 2026-05-25): when the input contains exactly one distinct radical `Sqrt[poly]` or `Power[poly, 1/q]` whose radicand is a polynomial in free symbols (e.g. `Sqrt[p+q]`, `Power[1+x^2, 1/3]`), `qa_cancel_with_poly_radical` substitutes the radical with a fresh symbol `S`, runs `Together`, reduces the numerator and denominator modulo `S^q - poly` via `PolynomialRemainder`, optionally rationalises via `PolynomialExtendedGCD` (when it shrinks the result), and substitutes back. Sample collapses: `Together[1/(x - Sqrt[p+q]) + 1/(x + Sqrt[p+q]), Extension -> Automatic]` returns `(2 x)/(-p - q + x^2)`; `Cancel[(x^2 - (p+q))/(x - Sqrt[p+q]), Extension -> Automatic]` returns `Sqrt[p+q] + x`. Multi-radical inputs (Cardano-style conjugate pairs) are rejected; the deeper limitation is documented in `docs/spec/changelog/2026-05-25.md`.

```mathematica
In[1]:= Together[a/b + c/d]
Out[1]= (b c + a d) / (b d)

In[2]:= Together[x^2/(x^2 - 1) + x/(x^2 - 1)]
Out[2]= x / (-1 + x)

In[3]:= Together[1/x + 1/(x + 1) + 1/(x + 2) + 1/(x + 3)]
Out[3]= (2 (3 + 11 x + 9 x^2 + 2 x^3)) / (x (1 + x) (2 + x) (3 + x))

In[4]:= Together[x^2/(x - y) - x y/(x - y)]
Out[4]= x

In[5]:= Together[y^(5/8)*(y^(19/8) - y^(73/24)/(y^(2/3) - 1/y^(1/3)))]
Out[5]= -y^3 / (-1 + y)

In[6]:= Together[1/(x - Sqrt[2]) + 1/(x + Sqrt[2]), Extension -> Sqrt[2]]
Out[6]= (2 x)/(-2 + x^2)
```

## Apart
Gives the partial fraction decomposition of a rational expression.
- `Apart[expr]`
- `Apart[expr, var]`
- `Apart[expr, var, Extension -> alpha]`

**Features**:
- `Protected`, `Listable`.
- Writes `expr` as a polynomial in `var` together with a sum of ratios of polynomials with minimal denominators.
- If `var` is not specified, intelligently selects the main polynomial variable natively.
- Implements exact undetermined coefficients algebraically leveraging row-reduced identity expansions over algebraic inputs avoiding recursive fractional losses natively.
- When `Together[expr]` produces a numerator or denominator that is not polynomial in the chosen variable (e.g. fractional-power inputs whose Together'd form is `y^(1/3)/(y - 1)`), the matrix-of-coefficients algorithm cannot apply; Apart returns the `Together` form unchanged rather than synthesising a spurious zero.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan) factors the denominator over `Q(alpha)` before partial-fraction decomposition runs, splitting reducible-over-extension factors (e.g. `x^2 - 2` into `(x - Sqrt[2])(x + Sqrt[2])` under `Extension -> Sqrt[2]`) and producing the corresponding linear-factor partial fractions. The pre-`Together` step is also extension-aware so any algebraic-number cancellations in numerator/denominator fire before splitting.

```mathematica
In[1]:= Apart[1/((1+x)(5+x))]
Out[1]= 1/(4 (1 + x)) - 1/(4 (5 + x))

In[2]:= Apart[(x^5-2)/((1+x+x^2)(2+x)(1-x))]
Out[2]= 2 - x + (-1 - x/3)/(1 + x + x^2) + 1/(9 (-1 + x)) - 34/(9 (2 + x))

In[3]:= Apart[(x+y)/((x+1)(y+1)(x-y)), x]
Out[3]= 2 y/((1 + y)^2 (x - y)) - (-1 + y)/((1 + x) (1 + y)^2)

In[4]:= Apart[1/(y^(2/3) - 1/y^(1/3))]
Out[4]= y^(1/3)/(-1 + y)

In[5]:= Apart[1/(x^2 - 2), x, Extension -> Sqrt[2]]
Out[5]= -1/2 1/(Sqrt[2] (Sqrt[2] + x)) + 1/4 Sqrt[2]/(-Sqrt[2] + x)
```

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

## GCD, LCM
- `GCD[n1, n2, ...]`: Greatest common divisor.
- `LCM[n1, n2, ...]`: Least common multiple.

Both fold through GMP when any argument is a bigint, so results that exceed `int64` (e.g. `LCM[20!, 10^100 + 3]`) are computed exactly rather than left symbolic.

## PowerMod
Gives modular exponentiations, inverses, and roots.
- `PowerMod[a, b, m]`: Gives $a^b \pmod m$.
- `PowerMod[a, -1, m]`: Finds the modular inverse of `a` modulo `m`.
- `PowerMod[a, 1/r, m]`: Finds a modular `r`-th root of `a`.

**Features**:
- `Protected`, `Listable`.
- Evaluates much more efficiently than `Mod[a^b, m]`.
- Integer-exponent path uses GMP `mpz_powm` / `mpz_invert`; `a`, `b`, `m`
  may all be arbitrary-precision bignums.
- Rational-exponent path `PowerMod[a, p/q, m]` solves the modular `q`-th
  root of `a^p` via Tonelli–Shanks, the coprime closed form, Hensel
  lifting, and CRT over the factorisation of `m`. The
  Adleman–Manders–Miller branch (`gcd(q, p-1) > 1` for `q > 2` over a
  large prime) is not yet supported.
- Returns unevaluated if the corresponding inverse or root does not exist
  (or cannot be computed by the implemented algorithms).
- Allows threading over lists natively.

```mathematica
In[1]:= PowerMod[2, 10, 3]
Out[1]= 1

In[2]:= PowerMod[3, -2, 7]
Out[2]= 4

In[3]:= PowerMod[3, 1/2, 2]
Out[3]= 1

In[4]:= PowerMod[2, {10, 11, 12, 13, 14}, 5]
Out[4]= {4, 3, 1, 2, 4}

In[5]:= PowerMod[100, 1/2, 17 * 19 * 23]
Out[5]= 10

In[6]:= PowerMod[2, 1/2, 10^18 + 9]
Out[6]= 742174169206529574
```

## PrimitiveRoot
Gives a primitive root of `n`, i.e. a generator of the multiplicative group
of integers modulo `n` that are coprime to `n`.
- `PrimitiveRoot[n]`: a primitive root of `n`.
- `PrimitiveRoot[n, k]`: the smallest primitive root of `n` greater than or
  equal to `k`.

**Features**:
- `Protected`, `Listable`.
- Returns unevaluated unless `n` is 2, 4, an odd prime power $p^k$, or
  twice an odd prime power $2 p^k$ (the moduli for which $(\mathbb{Z}/n\mathbb{Z})^*$
  is cyclic). For all other `n`, the call is left unevaluated.
- The 1-argument form returns a canonical primitive root: smallest for
  $n \in \{2, 4\}$ and odd prime powers; for $n = 2 p^k$ the formula
  $g$ if $g$ is odd else $g + p^k$ is applied, where $g$ is the smallest
  primitive root of $p^k$. This matches Mathematica's convention so that,
  e.g. `PrimitiveRoot[10] == 7` while `PrimitiveRoot[10, 1] == 3`.
- The 2-argument form walks forward from `k`; if `k > n - 1` the call is
  left unevaluated.
- All arithmetic uses GMP `mpz_t`, so machine integers, bignums, and
  symbolic bignum products like `Prime[1000000]^1000000` are handled
  uniformly. The prime-power detection iteratively strips prime exponents
  via `mpz_root`, which runs in $O(\omega(k))$ root extractions.
- Diagnostics:
  - `PrimitiveRoot::argt` if not called with 1 or 2 arguments.
  - `PrimitiveRoot::intg` if `n` (or the 2nd-arg `k` when numeric) is not
    an integer greater than 1.

```mathematica
In[1]:= PrimitiveRoot[9]
Out[1]= 2

In[2]:= PrimitiveRoot[10]
Out[2]= 7

In[3]:= PrimitiveRoot[10, 1]
Out[3]= 3

In[4]:= PrimitiveRoot[10, 4]
Out[4]= 7

In[5]:= PrimitiveRoot[{9, 7, 19}]
Out[5]= {2, 3, 2}

In[6]:= PrimitiveRoot[12]
Out[6]= PrimitiveRoot[12]
```

## PrimitiveRootList
Gives the sorted list of all primitive roots of `n` in the canonical
residues $\{1, \ldots, n-1\}$.
- `PrimitiveRootList[n]`

**Features**:
- `Protected`, `Listable`.
- Returns `{}` if `n` is not 2, 4, an odd prime power, or twice an odd
  prime power.
- Enumerates the $\varphi(\varphi(n))$ primitive roots as $g^i \bmod n$
  for $i \in \{1, \ldots, \varphi(n)\}$ with $\gcd(i, \varphi(n)) = 1$,
  where $g$ is the smallest primitive root of `n`. The list is sorted
  ascending.
- Falls back to unevaluated when $\varphi(n)$ exceeds `unsigned long`,
  since the enumeration cannot be represented.
- Non-integer numeric inputs (e.g. `11.0`, `11 + I`) flow through
  unevaluated with no diagnostic, matching Mathematica.
- Diagnostic: `PrimitiveRootList::argx` if not called with exactly 1
  argument.

```mathematica
In[1]:= PrimitiveRootList[9]
Out[1]= {2, 5}

In[2]:= PrimitiveRootList[19]
Out[2]= {2, 3, 10, 13, 14, 15}

In[3]:= PrimitiveRootList[12]
Out[3]= {}

In[4]:= Union[Table[PowerMod[2, i, 9], {i, 6}]]
Out[4]= {1, 2, 4, 5, 7, 8}
```

## MultiplicativeOrder
Gives the multiplicative order of `k` modulo `n` — the smallest positive
integer `m` such that $k^m \equiv 1 \pmod{n}$.
- `MultiplicativeOrder[k, n]`: order of `k` modulo `n`.
- `MultiplicativeOrder[k, n, {r1, r2, ...}]`: smallest positive `m` with
  $k^m \equiv r_i \pmod{n}$ for some `i` (a multi-target discrete log).

**Features**:
- `Protected`.
- All arithmetic uses GMP `mpz_t`, so `k`, `n`, and any `r_i` may be
  arbitrary-precision bignums.
- Negative `n` is treated as `|n|`; `k` (and each `r_i`) is reduced modulo
  `n` before the search, so negative or out-of-range inputs work
  transparently.
- Returns unevaluated when `gcd(k, n) != 1`, when `n` is zero, or when
  no power of `k` lands in the residue set (3-arg form).
- The 2-arg form factors `phi(n)` and successively strips prime factors
  from `phi(n)` whose corresponding exponent still maps `k` to 1 — so
  the work is dominated by factoring `phi(n)`, not by walking the orbit.
- The 3-arg form walks `k^m mod n` for `m = 1, ..., order(k, n)`. To
  guard against pathological group sizes, the call returns unevaluated
  if the order exceeds `10^8` or does not fit in an unsigned long.
- Non-integer numeric inputs (`Real`, `Complex`, `Rational`) and symbolic
  arguments flow through unevaluated with no diagnostic.
- Diagnostic: `MultiplicativeOrder::argt` when called with anything other
  than 2 or 3 arguments.

```mathematica
In[1]:= MultiplicativeOrder[5, 8]
Out[1]= 2

In[2]:= MultiplicativeOrder[5, 7]
Out[2]= 6

In[3]:= MultiplicativeOrder[-5, 7]
Out[3]= 3

In[4]:= MultiplicativeOrder[5, 7, {3, 11}]
Out[4]= 2

In[5]:= MultiplicativeOrder[10^10000, 7919]
Out[5]= 3959

In[6]:= Select[Range[43], MultiplicativeOrder[#, 43] == EulerPhi[43] &]
Out[6]= {3, 5, 12, 18, 19, 20, 26, 28, 29, 30, 33, 34}

In[7]:= MultiplicativeOrder[10, 22]
Out[7]= MultiplicativeOrder[10, 22]
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

## PrimeQ
- `PrimeQ[n]`: Returns `True` if `n` is a prime number, `False` otherwise.
- `PrimeQ[z]`: For a Gaussian integer `z = a + b I`, returns `True` if `z` is a Gaussian prime.
- `PrimeQ[n, GaussianIntegers -> True]`: Tests primality of `n` in the
  Gaussian integers `Z[i]`. A rational integer `n` is a Gaussian prime
  iff `|n|` is prime in `Z` AND `n ≡ 3 (mod 4)` (primes `≡ 1 (mod 4)`
  split as `(a+bi)(a-bi)`, and `2 = -i(1+i)^2` is reducible).

**Features**:
- `Listable`, `Protected`.
- Always returns `True` or `False`. For non-integer / non-Gaussian
  inputs (symbols, `Sqrt[2]`, `Exp[2 Pi I/3]`, strings, etc.) returns
  `False` — `*Q` predicates never remain symbolic.
- A Gaussian integer `a + b I` is a Gaussian prime if:
  - Both `a` and `b` are nonzero and `a^2 + b^2` is an ordinary prime, or
  - One of `a`, `b` is zero and the absolute value of the other is a prime congruent to 3 mod 4.

```mathematica
In[1]:= PrimeQ[7]
Out[1]= True

In[2]:= PrimeQ[1 + I]
Out[2]= True

In[3]:= PrimeQ[1 + 2 I]
Out[3]= True

In[4]:= PrimeQ[3 I]
Out[4]= True

In[5]:= PrimeQ[5 I]
Out[5]= False

In[6]:= PrimeQ[2 + 2 I]
Out[6]= False

In[7]:= PrimeQ[5, GaussianIntegers -> True]
Out[7]= False

In[8]:= PrimeQ[3, GaussianIntegers -> True]
Out[8]= True

In[9]:= PrimeQ[Exp[2 Pi I/3]]
Out[9]= False
```

## SquareFreeQ
- `SquareFreeQ[expr]`: Returns `True` if `expr` is a square-free polynomial or
  number, and `False` otherwise.
- `SquareFreeQ[expr, vars]`: Returns `True` if `expr` is square-free with
  respect to the variables `vars` (a `Symbol` or a `List` of `Symbol`s).
- `SquareFreeQ[..., GaussianIntegers -> True | False | Automatic]`: Tests
  square-freeness in `Z[i]` when `True`, in `Z` when `False`. `Automatic`
  (the default) enables Gaussian mode iff the input is a `Complex` literal
  with integer real and imaginary parts.

**Features**:
- `Protected`. Not `Listable` -- passing a list of inputs treats the list as
  the expression (`SquareFreeQ[{1, 2, 3}]` returns `False`).
- Always returns `True` or `False` on a structurally valid call. For inputs
  that are neither a recognised number nor a manifest polynomial -- reals,
  `Sqrt[2]`, `Sin[x]`, `Pi`, strings -- the result is `False`, never
  symbolic.
- An integer `n` is square-free iff `|n|` has no rational prime factor of
  multiplicity `>= 2`. `0` is not square-free; `+/-1` and `+/-p` (for any
  prime `p`) are.
- A rational `p/q` is square-free iff both `p` and `q` are square-free
  integers.
- A polynomial in `vars` is square-free iff for every variable `x_i` in
  `vars` that the polynomial actually depends on, `PolynomialGCD(p, dp/dx_i)`
  has degree `0` in `x_i`. Implementation routes the derivative through
  the `D` builtin and the gcd through `PolynomialGCD`.
- For `GaussianIntegers -> True` (or `Automatic` on a `Complex[Integer,
  Integer]` input), the test factors `N(z) = a^2 + b^2` over `Z` and
  classifies each rational prime by residue `mod 4`:
  - `p == 2`: the Gaussian prime above `2` is `1 + I`; its multiplicity in
    `z` equals the multiplicity of `2` in `N(z)`.
  - `p ≡ 3 (mod 4)`: `p` itself is a Gaussian prime with `N(p) = p^2`;
    its multiplicity in `z` is half the multiplicity of `p` in `N(z)`.
  - `p ≡ 1 (mod 4)`: `p` splits as `pi * conj(pi)`; ambiguous cases (the
    `e_in_norm == 2` slice) are resolved by computing `pi` via a
    Cornacchia search and stripping the multiplicity directly.
- The `Modulus -> p` option is parsed but only `Modulus -> 0` (the default
  no-modulus path) is wired in; any other value -- non-zero integer,
  non-integer, or symbolic -- emits `SquareFreeQ::modnotimpl` and leaves
  the call unevaluated until a polynomial sqfree-mod-`p` test is added.

Diagnostics:
- `SquareFreeQ[]` emits `SquareFreeQ::argb` and stays unevaluated.
- A non-`Rule` past the optional vars slot (e.g. `SquareFreeQ[1, 2, 3]`)
  emits `SquareFreeQ::nonopt` and stays unevaluated.
- `Modulus -> n` with `n != 0` emits `SquareFreeQ::modnotimpl` and stays
  unevaluated (e.g. `SquareFreeQ[x^2 + 1, Modulus -> 2]`).

```mathematica
In[1]:= SquareFreeQ[10]
Out[1]= True

In[2]:= SquareFreeQ[4]
Out[2]= False

In[3]:= SquareFreeQ[20]
Out[3]= False

In[4]:= SquareFreeQ[3 + 2 I]
Out[4]= True

In[5]:= SquareFreeQ[2, GaussianIntegers -> True]
Out[5]= False

In[6]:= SquareFreeQ[2/3]
Out[6]= True

In[7]:= SquareFreeQ[6 + 6 x + x^2]
Out[7]= True

In[8]:= SquareFreeQ[x^3 - x^2 y]
Out[8]= False

In[9]:= SquareFreeQ[x y^2, x]
Out[9]= True

In[10]:= SquareFreeQ[x y^2, y]
Out[10]= False

In[11]:= SquareFreeQ[10^70 + 3]
Out[11]= True
```

## IrreduciblePolynomialQ
- `IrreduciblePolynomialQ[poly]`: Returns `True` if `poly` is an irreducible
  polynomial over the rationals (treating any algebraic-number coefficient
  as an independent variable, per Mathematica's `Extension -> None` default).
- `IrreduciblePolynomialQ[poly, Extension -> alpha | {alpha_1, ...}]`: Tests
  irreducibility over the field extension `Q(alpha_1, alpha_2, ...)`.
- `IrreduciblePolynomialQ[poly, Extension -> Automatic]`: Extends `Q` by every
  algebraic-number coefficient appearing in `poly` (via `extension_autodetect`)
  and tests irreducibility over that compositum.
- `IrreduciblePolynomialQ[poly, Extension -> All]`: Tests absolute
  irreducibility over `C`.
- `IrreduciblePolynomialQ[poly, GaussianIntegers -> True]`: Tests irreducibility
  over `Q(i)`. If any coefficient of `poly` is a `Complex` literal, Gaussian
  mode is auto-enabled even without the option.

**Features**:
- `Listable`, `Protected`. The evaluator threads any `List` first argument
  element-wise; the builtin itself sees only scalar polynomial inputs.
- Always returns `True` or `False` on a structurally valid call. Constant
  numeric inputs (`0`, `1`, `5`, `2/3`), non-polynomial expressions (`Sin[x]`,
  `Pi`, `Sqrt[x]`), and the empty polynomial return `False`; constants are
  not irreducible polynomials.
- Algorithm: factor `poly` over the resolved field, then count non-constant
  factors with multiplicity. The polynomial is irreducible iff exactly one
  non-constant factor (mult `1`) appears in the factorisation. Numeric
  units (`Integer`, `Rational`, `Complex` leaves) and any sub-expression
  whose leaves are all free of the polynomial variables are treated as
  constants.
- Extension dispatch:
  - `Extension -> None` (default): `Factor[poly]` over `Q`. Atomic
    algebraic constants in `poly` (`Sqrt[int]`, `Power[int, p/q]`) are
    *frozen* to fresh placeholder symbols before factoring so Mathilda's
    multivariate Factor doesn't silently re-apply `Extension -> Automatic`
    -- matching Mathematica's "treat algebraic-number coefficients like
    independent variables" semantics.
  - `Extension -> alpha`: dispatches to the public `Factor[poly, Extension -> alpha]`.
  - `Extension -> {alpha_1, ...}`: list-of-generators tower, same routing.
  - `Extension -> Automatic`: `Factor[poly, Extension -> Automatic]` --
    Mathilda's own auto-detect.
  - `Extension -> All` (absolute irreducibility):
    * Univariate degree `1` -> `True`; degree `>= 2` -> `False` (splits
      over `C`).
    * Multivariate -- best-effort -- factors over `Q(i)` to catch
      conjugate-pair factorisations like `x^2 + y^2 = (x + i y)(x - i y)`.
- Gaussian mode (`GaussianIntegers -> True`, or auto-on for `Complex`-bearing
  input) calls `qa_factor_with_extension` directly on a `Complex[a, b] ->
  a + b I` lifted form of `poly` so the qa-factoring path sees the imaginary
  unit as a free symbol rather than an opaque `Complex` literal. The lifted
  form is intentionally un-evaluated to avoid Mathilda's
  `Times[b, I] -> Complex[0, b]` canonicalisation undoing the lift.

Multivariate inputs paired with an algebraic extension (`Extension -> α`,
`Extension -> All`, or `GaussianIntegers -> True`) take a separate path:
because the underlying `Factor[poly, Extension -> α]` only applies the
extension to single-variable polynomials, a multivariate poly that reaches
the no-extension fallback would otherwise look irreducible. To close that
gap, `IrreduciblePolynomialQ` runs a Hilbert-irreducibility specialisation
probe: pick the variable of maximum degree as the surviving univariate
indeterminate, substitute every other variable with each `c` in
`{2, 3, 5}`, factor the resulting univariate over the extension, and
flip the verdict to `False` only when *every* valid specialisation
produces `>= 2` non-constant factors. Hilbert's theorem says an
irreducible multivariate `p` stays irreducible under almost every integer
specialisation, so unanimous "reducible" probes are strong evidence; the
probe never downgrades a `False` from the cheap factor path (e.g. `x*y`).

Known limitations:
- The specialisation probe is necessary but not sufficient: a reducible
  multivariate polynomial whose factors degenerate at every probed `c`
  can slip through and report `True`. This is rare in practice and only
  affects polynomials with a small set of "bad" specialisations.
- `Modulus -> p` is not yet supported and is silently ignored.

Diagnostics:
- `IrreduciblePolynomialQ[]` emits `IrreduciblePolynomialQ::argx` and stays
  unevaluated.
- A trailing non-`Rule` past position `1` (e.g. `IrreduciblePolynomialQ[1, 2, 3]`,
  or an unknown option name like `IrreduciblePolynomialQ[x, Foo -> Bar]`)
  emits `IrreduciblePolynomialQ::nonopt` and stays unevaluated.

```mathematica
In[1]:= IrreduciblePolynomialQ[{x^2 - 1, x^2 - 2}]
Out[1]= {False, True}

In[2]:= IrreduciblePolynomialQ[{x^2 + 1, x^3 - 8}]
Out[2]= {True, False}

In[3]:= IrreduciblePolynomialQ[{x^4 - 4 y^2, x^4 - 2 y^2}]
Out[3]= {False, True}

In[4]:= IrreduciblePolynomialQ[x^2 + 2 I x - 1]
Out[4]= False

In[5]:= IrreduciblePolynomialQ[x^2 + 1, GaussianIntegers -> True]
Out[5]= False

In[6]:= IrreduciblePolynomialQ[x^2 + 2 Sqrt[2] x + 2]
Out[6]= True

In[7]:= IrreduciblePolynomialQ[x^2 + 2 Sqrt[2] x + 2, Extension -> Automatic]
Out[7]= False

In[8]:= IrreduciblePolynomialQ[{x^3 - 2, x^3 - 3}, Extension -> 2^(1/3)]
Out[8]= {False, True}

In[9]:= IrreduciblePolynomialQ[x^2 + 2, Extension -> {I, Sqrt[2]}]
Out[9]= False

In[10]:= IrreduciblePolynomialQ[{x^3 - 3, x^2 + 2 x y - 7}, Extension -> All]
Out[10]= {False, True}

In[11]:= IrreduciblePolynomialQ[{x^2 - 2 y^4, x^4 - 3 y^2}, Extension -> Sqrt[3]]
Out[11]= {True, False}

In[12]:= IrreduciblePolynomialQ[x^2 + y^2, GaussianIntegers -> True]
Out[12]= False

In[11]:= IrreduciblePolynomialQ[x^7 + 12 x y - 11, Extension -> All]
Out[11]= True
```

## PrimePi
- `PrimePi[x]`: Returns the number of primes less than or equal to `x`.

**Features**:
- `Listable`, `Protected`.
- Uses Meissel-Lehmer algorithm for efficient counting.

```mathematica
In[1]:= PrimePi[10]
Out[1]= 4

In[2]:= PrimePi[100]
Out[2]= 25

In[3]:= PrimePi[{10, 100}]
Out[3]= {4, 25}
```

## NextPrime
- `NextPrime[x]`: Smallest prime above `x`.
- `NextPrime[x, k]`: $k$-th next prime above `x` (or previous if $k$ is negative).

**Features**:
- `Protected`, `ReadProtected`.
- Supports negative $k$ for finding previous primes.
- Remains unevaluated if no such prime exists (e.g., `NextPrime[2, -1]`).

## FactorInteger
- `FactorInteger[n]`: Returns a list of prime factors and their exponents.
- `FactorInteger[n, k]`: Partial factorization, at most `k` distinct factors.
- `FactorInteger[n, Automatic]`: Pulls out easy factors using trial division.
- `FactorInteger[n, Method -> method]`: Factors `n` using a specific algorithmic method.

**Options for `Method`**:
- `"Automatic"`: (Default) Attempts factorization by sequentially executing Trial Division, Pollard's Rho, and ECM.
- `"TrialDivision"`: Extracts bounds matching the first 1000 computed primes and halts cleanly. Composite residues are preserved.
- `"PollardRho"`: Executes the Brent cycle-finding variant of Pollard's Rho algorithm targeting GMP bignums.
- `"BlakeRationalBaseDescent"`: Executes the Rational Base Descent algorithm against semiprimes `n = p q` where one factor is approximately `c (a/b)^k` for some unknown `c`, `k` and a small coprime integer pair `a > b`. With no explicit `"Base"`, auto-searches coprime partitions `a + b = j` for `j = 3, 4, ...` indefinitely until a factor is found (interruptible).
- `{"BlakeRationalBaseDescent", "Base" -> a/b}`: Same algorithm but pinned to a specific rational base `a/b`. Skips the auto-search.
- `"PollardP-1"`: Executes the Pollard $P-1$ algorithm, leveraging GMP and ECM natively.
- `"WilliamsP+1"`: Executes the Williams $P+1$ algorithm via the ECM library natively.
- `"Fermat"`: Explores factors symmetrically close to the square root boundary natively on large integers.
- `"CFRAC"`: Implements the continued fraction integer factorization method natively on GMP bignums.
- `"ShanksSquareForms"`: Implements Shanks's Square Forms Factorization (SQUFOF) natively on large integers. Halts explicitly if factors are not discovered within the loop constraints.
- `"ECM"`: Explicitly triggers Elliptic-Curve Method discovery natively.

**Features**:
- `Listable`, `Protected`.
- Supports negative integers (includes `{-1, 1}`).
- Supports rational numbers (denominator factors have negative exponents).

```mathematica
In[1]:= FactorInteger[12]
Out[1]= {{2, 2}, {3, 1}}

In[2]:= FactorInteger[-12]
Out[2]= {{-1, 1}, {2, 2}, {3, 1}}

In[3]:= FactorInteger[3/4]
Out[3]= {{2, -2}, {3, 1}}

In[4]:= FactorInteger[100, 1]
Out[4]= {{2, 2}, 25}
```

## EulerPhi
Gives the Euler totient function $\phi(n)$.
- `EulerPhi[n]`

**Features**:
- `Listable`, `Protected`.
- Counts the number of positive integers less than or equal to $n$ that are relatively prime to $n$.
- Returns 0 for $n = 0$, and handles negative integers via $\phi(-n) = \phi(n)$.
- Accepts arbitrary-precision integers (`BigInt`). Factorization runs in GMP
  through the same trial-division / Pollard-rho / ECM cascade used by
  `FactorInteger`, so inputs of cryptographic size are tractable.
- For a prime decomposition $n = \prod p_i^{k_i}$, computes
  $\phi(n) = n \prod (1 - 1/p_i)$ as $(n / \prod p_i) \prod (p_i - 1)$
  with exact integer arithmetic.

```mathematica
In[1]:= EulerPhi[10]
Out[1]= 4

In[2]:= EulerPhi[2^89 - 1]
Out[2]= 618970019642690137449562110
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

## Solve
Attempts to solve an equation or system of equations for one or more variables.
- `Solve[expr, vars]`: Solve `expr` for `vars` over the complex numbers (default).
- `Solve[expr, vars, dom]`: Solve over the domain `dom`. Supported: `Complexes` (default), `Reals`, `Integers`.

**Features**:
- `HoldAll`, `Protected`.
- Acts as a router that classifies its input and dispatches to a specialist:
  - Single equality, single variable -> `Solve`SolvePolynomialEquality` (below).
  - Single equality, single variable, polynomial specialist declines because
    the outermost dependence on `var` is an elementary invertible head ->
    inverse-function specialist (`src/solveinv.c`): peels `Log`, `Exp`,
    `Sin`/`Cos`/`Tan`/`Cot`/`Sec`/`Csc`, the hyperbolic counterparts, the
    inverse trig/hyperbolic forms, and `Power[g, n]` for integer `n >= 2`.
    Multi-branch heads introduce a fresh integer parameter `C[k]` and wrap
    each solution in `ConditionalExpression[..., Element[C[k], Integers]]`.
    Emits `Solve::ifun` on first use per call.
  - Single equality, single variable, both specialists above decline (because
    the equation carries `Sqrt[...]` / `x^(p/q)` / nested radicals) ->
    `Solve`SolveRadicalsEquality` (also below).
  - Multi-variable list, or `And`/`List` of equations -> `Solve`SolveLinearSystem`
    (also below).  The linear-system specialist accepts the same input shapes
    that the router uses to decide dispatch; it canonicalises each equation
    `lhs_i == rhs_i` to `lhs_i - rhs_i` and refuses (returns `NULL`) when the
    system is not affine in the variables, in which case the router leaves
    `Solve` unevaluated.
- Inequalities and multi-equation transcendental systems are reserved for
  future work and currently leave `Solve[...]` unevaluated.  When the
  inverse-function specialist's outermost peel succeeds but the inner
  equation is unsolvable and the peel was over `var` itself, Solve returns
  `{{var -> InverseFunction[head][rhs]}}` under `Solve::ifun`.
- **Approximate-number input**: if the equation contains any inexact numeric
  leaf (`Real` / MPFR), it is force-rationalised via the shared preprocessor
  in `src/common.c` before dispatch (so `1.5` becomes `3/2`, `N[Pi]` becomes
  a bit-exact rational, etc.), then the exact bindings produced by the
  specialist are numericalised on the way out -- same `inexact-in /
  inexact-out` contract `Integrate` and the exact-symbolic builtins
  (`Apart`, `Cancel`, `Together`, `Factor`, ...) follow.  The `vars`
  argument is never rationalised.  The preprocessor also tracks the
  *minimum* precision (in bits) across every inexact leaf and uses it
  both as the rationalisation tolerance and as the output precision, so a
  pure 30-digit-MPFR input flows back out at 30 digits, while a mixed
  Real + MPFR input drops to machine precision (the lower of the two)
  -- matching standard inexact-arithmetic semantics.
- Returns the solution set as a `List` of `List` of `Rule` pairs:
  - `{}` -- no solutions.
  - `{{}}` -- tautology (full-dimensional solution set).
  - `{{x -> v1}, {x -> v2}, ...}` -- one inner list per solution.  Multiplicity
    is preserved (repeated roots appear once per unit of multiplicity).
- **Rational-equality canonicalisation**: both sides are run through `Together`
  to combine into single fractions `N1/D1 == N2/D2`, then cross-multiplied to
  `N1*D2 - N2*D1 == 0` and `Collect`-ed in the solving variable before
  dispatch.  This routes equations like `a/x + b == 0` or `1/(x-1) == 2`
  through the polynomial specialist.  Any candidate root that provably zeroes
  one of the cleared denominators is dropped as extraneous (e.g.
  `Solve[x/(x-1) == 2/(x-1), x]` returns `{{x -> 2}}`, not `{{x -> 1}, {x -> 2}}`).
  Symbolic / undetermined denominator values are kept (parametric inputs like
  `Solve[a/x + b == 0, x]` return `{{x -> -a/b}}`).
- **Hidden-zero coefficient stripping**: after `Collect[Expand[...], var]` the
  per-degree coefficients are tested in turn with `PossibleZeroQ` (top down).
  Coefficients that test as zero but are not structurally zero -- e.g.
  `Sqrt[5 + 2 Sqrt[6]] - Sqrt[3] - Sqrt[2]`, recognised through the Stage-2
  numeric ladder -- are folded out and the polynomial is rebuilt at its true
  degree before the fast-path classifier sees it.  Without this pass the
  quadratic formula would divide by such a hidden-zero leading coefficient
  (`Solve[Sqrt[5 + 2 Sqrt[6]] x^2 - Sqrt[3] x^2 - Sqrt[2] x^2 - x - 1 == 0, x]`
  reduces to the linear `-x - 1 == 0` and returns `{{x -> -1}}`).  A
  hidden-zero constant is treated as a tautology
  (`Solve[Sqrt[5 + 2 Sqrt[6]] - Sqrt[3] - Sqrt[2] == 0, x]` returns `{{}}`).
- Per-degree handling for irreducible factors:
  - Degree 1 / 2: closed-form rules.
  - Quadratic in `Reals`: discriminant-aware.  Δ < 0 → no real roots;
    Δ = 0 → the double root is emitted *twice* (multiplicity preserved
    in step with the `Complexes` path); Δ > 0 → two distinct real
    roots.
  - Binomial `a*x^n + b == 0`: all n complex roots, or the real
    radical(s) in `Reals`.  Odd-`n` real branch: `(−b/a)^(1/n)` when
    `−b/a > 0`, `0` when `−b/a == 0`, and `−((b/a)^(1/n))` when
    `−b/a < 0` -- the last case is the *real* `n`-th root, not the
    principal complex one that `Power[base, 1/n]` produces by default.
    Even-`n`: ±r with `−b/a > 0`, `0` with `−b/a == 0`, `{}` with
    `−b/a < 0`.  Complex roots (no Reals constraint) are emitted as
    `r * (-1)^(2k/n)` for the principal radical `r = (-b/a)^(1/n)` and
    `k = 0..n-1`, then folded by `Power`'s rational-exponent canonicaliser
    so output matches Mathematica's standard form (e.g.
    `Solve[x^5 + 1 == 0, x]` returns
    `{{x -> (-1)^(1/5)}, {x -> (-1)^(3/5)}, {x -> -1},
       {x -> -(-1)^(2/5)}, {x -> -(-1)^(4/5)}}`).
  - n-quadratic `a*x^(2n) + b*x^n + c == 0`: substitution `u = x^n` followed by
    two binomial sub-solves; 2n radical roots regardless of `Cubics` / `Quartics`.
  - Degree 3: held `Root[Function[t, p[t]], k]` objects unless `Cubics -> True`.
  - Degree 4: held `Root[]` objects unless `Quartics -> True` (Ferrari is
    deferred -- with `Quartics -> True` the four roots are emitted as held
    `Root[]` for now).
  - Degree ≥ 5: held `Root[]` objects per irreducible factor.
- `Integers` domain is implemented as a post-pass over the `Reals` output:
  every candidate value is type-checked against `EXPR_INTEGER` /
  `EXPR_BIGINT` and dropped otherwise.  `Rational[p, q]`, irrational
  radicals (`Sqrt[2]`, `Power[2, 1/3]`, ...), held `Root[]` objects, and
  symbolic / parametric residues are *not* trusted to be integer-valued
  and are silently removed.  This means polynomials with one or more
  rational integer roots are returned correctly (`Solve[x^3 - 6 x^2 + 11
  x - 6 == 0, x, Integers]` -> `{{x -> 1}, {x -> 2}, {x -> 3}}` via
  factoring), but polynomials that only have irrational or symbolic
  integer roots return `{}`.  Higher-degree irreducibles default to
  `Root[]` form (`Cubics -> False`, `Quartics -> False`) and therefore
  yield `{}` under `Integers` unless the user opts into radical output.

**Options**:
- `Cubics -> False`: Emit cubic roots as held `Root[]` objects (default).
  `Cubics -> True` switches to closed-form Cardano radicals.
- `Quartics -> False`: Emit quartic roots as held `Root[]` objects (default).
  Reserved: Ferrari closed form is deferred.
- `InverseFunctions -> Automatic`: Enables the inverse-function specialist
  (default).  Set to `False` to disable the specialist; equations that can
  only be solved through inversion then return unevaluated.
- `GeneratedParameters -> C`: Head used by the inverse-function specialist
  when minting fresh integer-parameter symbols `C[1], C[2], ...`.  Only the
  bare-symbol form is honoured; the `Function` form is reserved.
- `VerifySolutions -> Automatic`: Reserved.

```mathematica
In[1]:= Solve[2 x + 3 == 0, x]
Out[1]= {{x -> -3/2}}

In[2]:= Solve[x^2 - 5 x + 6 == 0, x]
Out[2]= {{x -> 2}, {x -> 3}}

In[3]:= Solve[x^2 + 1 == 0, x]
Out[3]= {{x -> -I}, {x -> I}}

In[4]:= Solve[x^2 + 1 == 0, x, Reals]
Out[4]= {}

In[5]:= Solve[(x-1)^2 == 0, x]
Out[5]= {{x -> 1}, {x -> 1}}

In[6]:= Solve[x^4 - 5 x^2 + 4 == 0, x]
Out[6]= {{x -> 1}, {x -> -1}, {x -> 2}, {x -> -2}}

In[7]:= Solve[x^3 + x + 1 == 0, x]
Out[7]= {{x -> Root[Function[1 + #1 + #1^3], 1]}, ...}

In[8]:= Solve[Sin[x] == 0, x]
Out[8]= Solve[Sin[x] == 0, x]

In[9]:= Solve[a/x + b == 0, x]
Out[9]= {{x -> -a/b}}

In[10]:= Solve[1/(x-1) == 2, x]
Out[10]= {{x -> 3/2}}

In[11]:= Solve[x/(x-1) == 2/(x-1), x]
Out[11]= {{x -> 2}}             (* x = 1 dropped as extraneous *)

In[12]:= Solve[x^2 - 5 x + 6 == 0, x, Integers]
Out[12]= {{x -> 2}, {x -> 3}}

In[13]:= Solve[x^2 - 2 == 0, x, Integers]
Out[13]= {}                     (* Sqrt[2] is not an Integer *)

In[14]:= Solve[1.5 x + 3 == 0, x]
Out[14]= {{x -> -2.0}}          (* approximate-in / approximate-out *)

In[15]:= Solve[{1.5 x + y == 4.5, x - y == 0.5}, {x, y}]
Out[15]= {{x -> 2.0, y -> 1.5}}

In[16]:= Solve[N[Pi, 50] x == 1, x]
Out[16]= {{x -> 0.31830988618379067...}}   (* 50-digit MPFR result *)

In[14]:= Solve[x^3 + 1 == 0, x, Reals]
Out[14]= {{x -> -1}}            (* real cube root, not (-1)^(1/3) *)

In[15]:= Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x, Integers]
Out[15]= {{x -> 1}, {x -> 2}, {x -> 3}}

In[16]:= Solve[3 x + 2 y == 11 && x + y == 12, {x, y}]
Out[16]= {{x -> -13, y -> 25}}

In[17]:= Solve[a x + c == 1 && b x - d y == 2, {x, y}]
Out[17]= {{x -> (1 - c)/a, y -> (-2 a + b - b c)/(a d)}}

In[18]:= Solve[3 x + 2 y == 11, {x, y}]
Solve::svars: Equations may not give solutions for all "solve" variables.
Out[18]= {{y -> 11/2 - 3 x/2}}     (* x is free; only y has a rule *)

In[19]:= Solve[3 x + 2 y == 11 && x + y == 12 && 3 x + y == 32, {x, y}]
Out[19]= {}                        (* over-determined, inconsistent *)
```

## Solve`SolveLinearSystem
The linear-system specialist invoked by `Solve` for multi-variable inputs
(`And` / `List` of equations, or a single equation paired with a multi-symbol
variable list).  Reachable directly via its context-qualified name when the
caller has already classified its input.
- `Solve`SolveLinearSystem[eqns, vars]`
- `Solve`SolveLinearSystem[eqns, vars, dom]`

**Features**:
- `Protected`.
- `eqns` may be a single `Equal[lhs, rhs]`, `And[Equal[...], ...]`, or
  `List[Equal[...], ...]`.  `vars` must be a `List` of distinct symbols.
- Each equation is canonicalised to `lhs - rhs` and `Expand`-ed, then
  asserted affine in `vars`: coefficient of each `var` must be free of the
  variables, and the residual after subtracting `sum_j coeff_j * vars[j]`
  must also be free of the variables.  Non-affine systems return `NULL`
  (caller leaves `Solve` unevaluated).
- The m x (n+1) augmented matrix is built with variable columns in
  **reversed order** (M[i][0] is the coefficient of `vars[n-1]`).  This is
  what produces Mathematica's `Solve::svars` convention for under-determined
  systems: left-to-right Gauss--Jordan then naturally pivots on the
  right-most variable first, leaving left-most variables free.
- Gauss--Jordan elimination with symbolic-pivot selection: among non-zero
  candidates in the current column, prefer a concretely-non-zero entry
  (`Integer`, `Rational`, `Real`) over a symbolic one.  A column whose
  entries all simplify to zero (via `Cancel[Together[.]]`) becomes a free
  variable.  After reduction, any zero row whose augmented column is
  non-zero is detected as an inconsistency.
- Output shape:
  - Unique solution: `{{v1 -> e1, v2 -> e2, ...}}` (rules in input order).
  - Inconsistent system: `{}`.
  - Under-determined system: `{{pivot_vars -> exprs in free vars}}`; emits
    `Solve::svars`; free variables produce no rule.
  - Empty equation list (`Solve[True, vars]`): `{{}}` (tautology).
- Domain filtering (post-pass):
  - `Integers`: every emitted rule's RHS must be a concrete `EXPR_INTEGER`;
    otherwise the entire solution is dropped (`{}`).
  - `Reals`: any RHS that syntactically contains a `Complex[_, _]` head is
    treated as non-real and the whole solution is dropped.
  - `Complexes` (default): no filter.

```mathematica
In[1]:= Solve`SolveLinearSystem[{x + y == 3, x - y == 1}, {x, y}]
Out[1]= {{x -> 2, y -> 1}}

In[2]:= Solve`SolveLinearSystem[{x + y == 0}, {x, y}]
Solve::svars: Equations may not give solutions for all "solve" variables.
Out[2]= {{y -> -x}}                (* x free *)

In[3]:= Solve`SolveLinearSystem[
            {x + y + z == 6, 2 x - y + z == 3, x + 2 y - z == 2},
            {x, y, z}]
Out[3]= {{x -> 1, y -> 2, z -> 3}}
```

## Solve`SolvePolynomialEquality
The polynomial-equality specialist invoked by `Solve`.  Reachable directly via
its context-qualified name when the caller has already classified its input.
- `Solve`SolvePolynomialEquality[lhs == rhs, var]`
- `Solve`SolvePolynomialEquality[lhs == rhs, var, dom]`

**Features**:
- `Protected`.
- Same algorithm and output shape as `Solve` for single polynomial equalities
  in one variable.  Does not parse options; the caller supplies them through
  the C-level entry point.

## Solve`SolveInverseFunctions
The inverse-function specialist invoked by `Solve` when the outermost
dependence on `var` is an elementary invertible head.  Reachable directly
via its context-qualified name when the caller has already classified its
input.
- `Solve`SolveInverseFunctions[lhs == rhs, var]`
- `Solve`SolveInverseFunctions[lhs == rhs, var, dom]`

**Features**:
- `Protected`.
- Recognised heads: `Log`, `Exp`, `Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc`,
  the hyperbolic counterparts (`Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`,
  `Csch`), their inverses (`ArcSin`, `ArcCos`, ..., `ArcCsch`), and
  `Power[g, n]` for integer `n >= 2`.  Also recognises `Power[E, g(x)]`
  as the canonical form of `Exp[g(x)]`.
- Additive-shift isolation pre-pass: equations of the form
  `c * head[g(x)] + free_of_var == 0` are reduced to
  `head[g(x)] == new_rhs` before head dispatch.
- Multi-branch heads introduce a fresh integer parameter `C[k]` (head
  controlled by the parent `Solve`'s `GeneratedParameters` option) and
  wrap each branch in `ConditionalExpression[..., Element[C[k], Integers]]`.
- Inverse heads (`ArcSin`, `ArcCos`, `ArcTan`) use a vertical-strip
  predicate on `Re[a]`/`Im[a]` matching Mathematica's principal-branch
  domain.
- Inner equations are solved by hand-off to
  `Solve`SolvePolynomialEquality` -> `Solve`SolveInverseFunctions`
  (depth-capped at 8) -> `Solve`SolveRadicalsEquality`.  The recursion
  bypasses `Solve` itself, so the parent's inexact-rationalisation pre-
  pass runs only once.
- Emits `Solve::ifun` to stderr on first multi-branch peel per call and
  on the `InverseFunction[head][rhs]` fallback.
- Does not parse options; the caller supplies them through the C-level
  entry point `solveinv_solve_inverse_equality`.  When called via the
  qualified builtin, defaults `InverseFunctions -> True` and
  `GeneratedParameters -> C` are used.

```mathematica
In[1]:= Solve`SolveInverseFunctions[Sin[x] == a, x]
Out[1]= {{x -> ConditionalExpression[Pi - ArcSin[a] + 2 Pi C[1],
                                     Element[C[1], Integers]]},
         {x -> ConditionalExpression[ArcSin[a] + 2 Pi C[1],
                                     Element[C[1], Integers]]}}

In[2]:= Solve`SolveInverseFunctions[Log[x^2 + 1] + 1 == 0, x]
Out[2]= {{x -> -I Sqrt[1 - 1/E]}, {x -> I Sqrt[1 - 1/E]}}
```

## Solve`SolveRadicalsEquality
The radicals-equation specialist invoked by `Solve` when the polynomial
specialist declines (because the equation contains `Sqrt[...]`, fractional
powers `x^(p/q)`, or nested radicals).  Reachable directly via its
context-qualified name when the caller has already classified its input.
- `Solve`SolveRadicalsEquality[lhs == rhs, var]`
- `Solve`SolveRadicalsEquality[lhs == rhs, var, dom]`

**Features**:
- `Protected`.
- Algorithm:
  1. Canonicalise the equation: compute `e = Numerator[Together[lhs - rhs]]`.
     This combines rational-radical inputs (e.g.
     `(x + Sqrt[x])/Sqrt[x] + Sqrt[x]/(x + Sqrt[x]) == 4`) into a single
     polynomial-style residual in `var` and the radicals it contains.
  2. Iteratively locate radical atoms `Power[base, p/q]` (q > 1) anywhere in
     the working system (main equation + accumulated side equations).  For
     each distinct base `g_i`, introduce a fresh generator `u_i` so that
     `u_i = g_i^(1/L_i)`, where `L_i` is the LCM of denominators of *all*
     exponents of `g_i` (so `Sqrt[x]` and `x^(1/4)` share a single
     generator `u = x^(1/4)` with `L = 4`).  Replace every `g_i^(p/q)`
     by `u_i^(p*L_i/q)` in every equation, and append the side equation
     `u_i^L_i - g_i == 0`.  Nested radicals are picked up automatically
     -- a fresh atom inside a previously substituted base becomes its
     own generator on the next iteration.
  3. Eliminate `u_1, u_2, ...` from the main equation by chained
     `Resultant_{u_i}(main, side_eq_i, u_i)` in introduction order, so
     each side equation contributes exactly one fresh generator and the
     end-result is a polynomial in `var` alone.
  4. Hand the eliminated polynomial to `Solve`SolvePolynomialEquality`.
  5. Verify every candidate by back-substitution into the *original*
     equation.  The residual is first evaluated numerically with `N[]`
     and rejected when its magnitude exceeds `1e-9`; only when the
     numerical pass cannot decide -- free parameters, removable-
     singularity `Indeterminate`, etc. -- does the verifier fall back
     to a symbolic `Simplify` pass to catch structural zeros.  This
     ordering matters for candidates with algebraic coefficients
     (e.g. `Sqrt[2]` in the elimination): `Simplify` on the back-
     substituted residual can run for seconds per candidate, while
     `N[]` evaluates the same residual in microseconds.  Candidates
     whose residual still depends on free parameters (and so cannot
     be decided either way) are kept and trigger `Solve::nongen`,
     matching Mathematica's convention.
- Output shape matches `Solve`SolvePolynomialEquality`: a `List` of
  singleton-rule `List`s, plus the empty `List[]` when no candidate
  survives verification.  The `dom` argument flows through to the
  polynomial specialist (so `Reals` filters the candidate polynomial
  via the same per-degree discriminant tests as the polynomial path).
- The substitution introduces fresh generator symbols whose names follow
  the template `$radu<n>$`.  They are local to the call -- they never
  appear in the returned solution list (the resultant elimination
  removes every generator).
- The verifier accepts `Root[poly, k]` candidates without further
  checks: the polynomial elimination is exact, and `Root[]` objects
  describe the unique algebraic root of an irreducible factor that
  is not amenable to back-substitution.  Reflects Mathematica's
  policy of keeping `Root[]`-form solutions when they cannot be
  further simplified.
- The substitution-then-elimination strategy is "complete up to
  verification": every actual solution survives if it is closed-form,
  while extraneous roots from cross-multiplication (Together) or
  L-th-root branching are filtered out at the verifier.

```mathematica
In[1]:= Solve[Sqrt[x] + 3 x == 5, x]
Out[1]= {{x -> 1/18 (31 - Sqrt[61])}}

In[2]:= Solve[Sqrt[x] + 3 == 5, x]
Out[2]= {{x -> 4}}

In[3]:= Solve[x - 8 Sqrt[x] + 15 == 0, x]
Out[3]= {{x -> 9}, {x -> 25}}

In[4]:= Solve[Sqrt[x] + 3 x^(1/4) == 5, x]
Out[4]= {{x -> 1/2 (311 - 57 Sqrt[29])}}

In[5]:= Solve[(x + Sqrt[x])/Sqrt[x] + Sqrt[x]/(x + Sqrt[x]) == 4, x]
Out[5]= {{x -> 2 (2 + Sqrt[3])}}

In[6]:= Solve[Sqrt[x + 1] + Sqrt[x - 1] == 3, x]
Out[6]= {{x -> 85/36}}

In[7]:= Solve[Sqrt[x + 5] + Sqrt[x] == -1, x]
Out[7]= {}

In[8]:= Solve[x + Sqrt[x - 1] == 1, x]
Out[8]= {{x -> 1}}

In[9]:= Solve[Sqrt[a x + c] + 3 x == 5, x]
Solve::nongen: There may be values of the parameters for which some or
              all solutions are not valid.
Out[9]= {{x -> 1/18 (30 + a - Sqrt[60 a + a^2 + 36 c])},
         {x -> 1/18 (30 + a + Sqrt[60 a + a^2 + 36 c])}}
```

## Cubics
Option for `Solve` that controls whether cubic equations are solved via
explicit radical formulas.
- `Cubics -> False` (default): emit held `Root[]` objects.
- `Cubics -> True`: emit closed-form Cardano radicals.

## Quartics
Option for `Solve` that controls whether quartic equations are solved via
explicit radical formulas.
- `Quartics -> False` (default): emit held `Root[]` objects.
- `Quartics -> True`: reserved -- Ferrari is deferred; held `Root[]` objects
  are still emitted in the current implementation.

