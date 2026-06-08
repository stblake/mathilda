# Algebra

Algebraic manipulation of rational and polynomial expressions and equation solving: rational-expression normal forms (`Numerator`, `Denominator`, `Cancel`, `Together`, `Apart`), resultant-based tools, minimal and irreducible-polynomial predicates, and the equation solver (`Solve`, `Root`, `SolveAlways`, `ToRadicals`, and the cubic/quartic closed-form helpers).

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

## Subresultants

- `Subresultants[poly1, poly2, var]`: Gives the list of **principal
  subresultant coefficients** (PSCs) of `poly1` and `poly2` treated as
  polynomials in `var`.
- The list has length `Min[Exponent[poly1, var], Exponent[poly2, var]] + 1`.
- The first element equals `Resultant[poly1, poly2, var]`; the first `k`
  entries vanish exactly when the polynomials have `k` common roots
  (multiplicity counted), so the index of the first non-zero entry is the
  degree of `PolynomialGCD[poly1, poly2, var]`.
- Computed efficiently by a subresultant polynomial-remainder sequence (the
  same Bronstein γ/β/δ recurrence as `Resultant`), reading off the leading
  coefficient `s_p = lc(R_p)^{δ_p} / s_{p-1}^{δ_p - 1}` at each chain degree
  (the cumulative `s_{p-1}` factor is what makes degree gaps correct even when
  the preceding coefficient is not a unit); defective (degree-gap) indices come
  out as `0`. Inputs with algebraic-number coefficients (e.g. `Sqrt[2]`) fall
  back to the Sylvester-minor determinant definition. `var` must be a symbol;
  PSCs are polynomials in the coefficients of the inputs.

```
In[1]:= Subresultants[2x^7 + 3x^3 - 7x + 1, 3x^5 - 17x + 21, x]
Out[1]= {273612691817, 68946901, 1299537, 16641, 0, 9}

In[2]:= Subresultants[(x - 1)(x - 2)(x - 3), (x - 1)(x - 4)(x - 5), x]
Out[2]= {0, 12, -4, 1}

In[3]:= Subresultants[(x - 1)^5 (x - 2)(x - 3), (x - 1)^4 (x - 4)(x - 5), x]
Out[3]= {0, 0, 0, 0, 144, 18, 1}

In[4]:= Subresultants[a x^3 + b x^2 + c x + d, x^3 - 5 b x - 7 a, x]
Out[4]= {-343 a^6 + ... - d^3, -7 a^2 b + 25 a^2 b^2 - 5 b^3 + 10 a b c + c^2 - b d, -b, 1}

In[5]:= Length[Subresultants[x^50 + a, x^20 + b, x]]
Out[5]= 21
```

## SubresultantPolynomials

- `SubresultantPolynomials[poly1, poly2, var]`: Gives the list of
  **subresultant polynomials** `{S_0, ..., S_m}` of `poly1` and `poly2` treated
  as polynomials in `var`, where `m = Exponent[poly2, var]`.
- Requires `Exponent[poly1, var] >= Exponent[poly2, var]` and **exact**
  coefficients; otherwise emits `SubresultantPolynomials::npolys` and stays
  unevaluated.
- The list has length `m + 1`. Its first element is
  `Resultant[poly1, poly2, var]`; its last element is `poly2` (up to a leading
  power of `lc[poly2]` when `Exponent[poly1, var] > Exponent[poly2, var] + 1`).
- The degree of `S_j` is at most `j`, and the coefficient of `var^j` in `S_j`
  is the `j`-th principal subresultant coefficient, i.e.
  `Subresultants[poly1, poly2, var][[j+1]]`. Defective (degree-gap) entries are
  lower-degree polynomials or `0`.
- Computed efficiently by the subresultant polynomial-remainder sequence (the
  same Bronstein γ/β/δ recurrence as `Subresultants`): each chain member `R_p`
  is the defective subresultant of index `deg(R_{p-1}) - 1`, and the regular
  subresultant of degree `deg(R_p)` is `R_p` rescaled by
  `(lc(R_p) / s_prev)^{δ_p - 1}`. Algebraic-number coefficients fall back to the
  determinant-polynomial definition.

```
In[1]:= SubresultantPolynomials[(x - 1)^2 (x - 2) (x - 3), (x - 1) (x - 4)^2, x]
Out[1]= {0, -36 + 36 x, 38 - 49 x + 11 x^2, -16 + 24 x - 9 x^2 + x^3}

In[2]:= SubresultantPolynomials[a x^3 + b x^2 + c x + d, 3 a x^2 + b x + c, x]
Out[2]= {4 a^2 c^3 + 2 a b^3 d - 18 a^2 b c d + 27 a^3 d^2,
         -2 a b c + 9 a^2 d - 2 a b^2 x + 6 a^2 c x, c + b x + 3 a x^2}

In[3]:= SubresultantPolynomials[2 x^7 + 3 x^3 + 5 x - 1, 7 x^6 + 8 x - 9, x]
Out[3]= {-183782157189, -761749829 + 3208696817 x,
         -3143546 + 11222638 x + 3838135 x^2,
         -21609 + 163611 x - 49392 x^2 + 64827 x^3, 0,
         -49 + 371 x - 112 x^2 + 147 x^3, -9 + 8 x + 7 x^6}

In[4]:= First[%] - Resultant[2 x^7 + 3 x^3 + 5 x - 1, 7 x^6 + 8 x - 9, x]
Out[4]= 0
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

## MinimalPolynomial

- `MinimalPolynomial[s, x]`: Returns the lowest-degree polynomial in `x` with
  integer coefficients, positive leading coefficient, and content `1` (GCD of
  coefficients equal to `1`) having the algebraic number `s` as a root.
- `MinimalPolynomial[s]`: Returns the minimal polynomial as a pure function
  (e.g. `1 - 10 #1^2 + #1^4 &`).
- `MinimalPolynomial[s, x, Extension -> a]`: Returns the characteristic
  polynomial of `s` in `Q(a)` over `Q(a)` — equal to `m_s(x)^([Q(a):Q]/[Q(s):Q])`
  when `s` lies in `Q(a)`; otherwise the call is left unevaluated.

**Features**:
- `Listable`, `Protected`. A `List` first argument threads element-wise.
- `s` may be built from integers and rationals, radicals (`Sqrt`,
  `Power[_, p/q]`), the imaginary unit, roots of unity (`Power[E, I Pi r]`),
  `Root[f &, k]` objects, and the field operations `Plus`/`Times`/`Power`.
- Non-algebraic input (e.g. `Pi`, `Log[2]`) is left unevaluated.

**Algorithm**: each algebraic atom of `s` is replaced by a fresh auxiliary
variable carrying a polynomial defining relation (negative powers become
reciprocal variables `D w - 1`, so every relation stays polynomial). The
auxiliaries are eliminated from `x - s` by repeated `Resultant`, the result is
made primitive over `Z` and factored, and the unique irreducible factor that
vanishes at `s` (chosen by a high-precision numeric test) is returned. The
`Extension` form uses the tower law on the degrees produced by the same core,
with membership `s ∈ Q(a)` verified through the primitive-element degree
`[Q(a, s) : Q] = [Q(a) : Q]`.

```
In[1]:= MinimalPolynomial[Sqrt[2] + Sqrt[3], x]
Out[1]= 1 - 10 x^2 + x^4

In[2]:= MinimalPolynomial[(1 + I)/Sqrt[2], x]
Out[2]= 1 + x^4

In[3]:= MinimalPolynomial[Root[2 #1^3 - 2 #1 + 7 &, 1] + 17, x]
Out[3]= -9785 + 1732 x - 102 x^2 + 2 x^3

In[4]:= MinimalPolynomial[Sqrt[2], x, Extension -> E^(I Pi/4)]
Out[4]= 4 - 4 x^2 + x^4
```

## Solve

Attempts to solve an equation or system of equations for one or more variables.
- `Solve[expr, vars]`: Solve `expr` for `vars` over the complex numbers (default).
- `Solve[expr, vars, dom]`: Solve over the domain `dom`. Supported: `Complexes` (default), `Reals`, `Integers`.

**Features**:
- `Protected`.  Matches Mathematica's attribute set -- arguments are
  evaluated by the evaluator before reaching the router.  When the
  second argument has been substituted to a numeric atom (typically
  because an OwnValue like `x = 5` was previously set, or the user
  literally passed `Solve[..., 5]`), the router emits `Solve::ivar`
  and returns unevaluated.
- **Generalised (compound) variables.**  `vars` may contain any
  non-numeric expression, not only symbols: `Solve[lhs == rhs, Dt[y]]`,
  `Solve[f[a] + b == c, f[a]]`, `Solve[a x^2 + b == 0, x^2]`, and
  multi-var forms like `Solve[{...}, {Dt[x], Dt[y]}]` are all
  accepted.  The router substitutes each non-symbol entry with a
  fresh internal symbol (`Solve$var$N`), runs the standard dispatch,
  and reverses the substitution on the result so the user sees Rule
  LHSes like `Dt[y] -> ...` directly.  The substitution is purely
  structural (literal `expr_eq`); polynomial identifications like
  `x^4 == (x^2)^2` are not yet recognised, so `Solve[x^4 - 1 == 0,
  x^2]` returns the substituted form rather than `{{x^2 -> 1},
  {x^2 -> -1}}`.
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

## Root


Held symbolic representation of an indexed root of a univariate polynomial.

- `Root[Function[t, p[t]], k]` — the `k`-th root of `p` (1-indexed).

**Canonical index `k`** (matches Mathematica):

1. **Real roots first**, ordered ascending by value.
2. **Complex roots** afterwards, ordered by `Re` ascending; ties broken by
   `|Im|` ascending; within a conjugate pair the negative-`Im` member comes
   first.

This is the convention used by both `Solve`'s emission and `N[Root[..]]`'s
numericalization, so `Root[f, 1]` always refers to the same root regardless
of how it was produced.

**Numericalization** — `N[Root[f, k]]` and `N[Root[f, k], prec]`:

The pipeline is companion-matrix QR → Sturm certificate → canonical sort →
Newton refinement → basin verification. Both real and complex roots are
returned as MPFR values (`Complex[MPFR, MPFR]` for complex). Failure modes:

- `Root::nonint` — polynomial has non-integer coefficients (deferred case).
- `Root::indx`   — `k` is outside `1..deg(p)`.
- `Root::conv`   — QR or Newton did not converge after one precision retry.

Examples:
```
In[1]:= N[Root[Function[#^3 - 2 # - 5], 1], 30]
Out[1]= 2.094551481542326591482386540580

In[2]:= N[Root[Function[#^3 + # + 1], 1], 20]    (* real root first *)
Out[2]= -0.68232780382801932737

In[3]:= N[Root[Function[#^3 + # + 1], 2], 20]    (* conj pair: -Im first *)
Out[3]= 0.34116390191400966368 - 1.1615414252683233453 I

In[4]:= N[Root[Function[#^3 + # + 1], 3], 20]
Out[4]= 0.34116390191400966368 + 1.1615414252683233453 I
```

## ToRadicals


Convert held `Root[Function[poly], k]` objects in an expression into
closed-form radical expressions.

- `ToRadicals[expr]`

**Features**:

- `Protected`.
- Closed-form radicals are always returned when the polynomial has degree
  at most four — linear (trivial), quadratic (`Sqrt`), cubic (Cardano), and
  quartic (Ferrari via the depressed quartic + resolvent cubic).
- Binomial Root objects `Root[Function[a #^n + b], k]` are reduced to
  radicals for any degree `n`, using the principal `n`-th root multiplied
  by `(-1)^(2 (k-1) / n)`.
- Other Root objects of degree ≥ 5 are returned unchanged — the system
  makes no attempt at decomposition or solvable-Galois detection (cf.
  Mathematica's note "ToRadicals cannot find them").
- The k-th radical root is selected to agree with `N[Root[poly, k]]`'s
  canonical ordering (real-first ascending, complex by `Re` / `|Im|` /
  negative-`Im` first) — each formula's natural emission order is
  numerically matched against `root_numericalize` at machine precision.
  When the polynomial carries parametric coefficients (no numericalisation
  possible), the natural per-formula index `k - 1` is used and the result
  is allowed to disagree with `expr` for some parameter values, matching
  Mathematica's `nongen` behaviour.
- Walks its argument recursively, so `Root[..]` nodes inside `List`,
  `Equal`, `Less`, `Greater`, `And`, `Or`, `Not`, `Implies`, ... thread
  automatically — every `Root` anywhere in the tree is processed
  independently and the surrounding structure is preserved.
- Idempotent: `ToRadicals[ToRadicals[expr]] === ToRadicals[expr]`, since a
  successful conversion produces an expression free of `Root[..]` nodes.

```mathematica
In[1]:= ToRadicals[Root[Function[#^2 + 3 # + 5], 1]]
Out[1]= 1/2 (-3 - I Sqrt[11])

In[2]:= ToRadicals[Root[Function[#^2 + 3 # + 5], 2]]
Out[2]= 1/2 (-3 + I Sqrt[11])

In[3]:= ToRadicals[Root[Function[#^5 - 2], 3]]
Out[3]= (-1)^(4/5) 2^(1/5)

In[4]:= With[{r = ToRadicals[Root[Function[#^4 + 3 #^3 - 5 #^2 - 7 # + 9], 1]]},
              Chop[N[r^4 + 3 r^3 - 5 r^2 - 7 r + 9, 30]]]
Out[4]= 0

In[5]:= ToRadicals[Root[Function[#^5 - # - 1], 1]]      (* non-binomial deg 5 *)
Out[5]= Root[#1^5 - #1 - 1 &, 1]

In[6]:= ToRadicals[Root[Function[#^2 - 2], 2] < 3]      (* threading *)
Out[6]= True
```

## SolveAlways

Finds values of parameters (the symbols appearing in `eqns` but **not** in
`vars`) that make every equation in `eqns` hold for every value of `vars`.
The reduction is: each `lhs == rhs` is rewritten as the polynomial
`lhs - rhs`; `CoefficientList[lhs - rhs, vars]` exposes every coefficient
as a polynomial in the remaining symbols; every such coefficient must
vanish; the resulting system is fed to `Solve` with the parameters as
unknowns.
- `SolveAlways[eqns, vars]`

**Scope**:
- `eqns` may be `Equal[lhs, rhs]`, a `List[Equal[...], ...]`, or an
  `And[Equal[...], ...]`.  `True`/`False` sentinels arising from the
  evaluator's pre-pass on `==` are folded into the trivial answers.
- `vars` is a symbol or a `List` of symbols.
- The empty-parameter case (every symbol in `eqns` already appears in
  `vars`) returns `{}` — Mathematica's convention that there are no
  parameter values to report regardless of whether the polynomial is
  identically zero.
- `Unequal` (`!=`), `Or`-combinations of equations, radicals
  (`Sqrt[a x] == ...`), and `Series`-with-`O[x]^n` stripping are **not**
  handled in this version; those inputs will produce a `Solve`-level
  result on whatever coefficient system `CoefficientList` produces, which
  may not be the SolveAlways-correct answer.

**Diagnostics**:
- `SolveAlways::argt` — wrong number of arguments.
- `SolveAlways::eqf` — `eqns` contained a non-`Equal` element.
- `SolveAlways::ivar` — `vars` was not a symbol or non-empty list of
  symbols.

**Examples**:

```
In[1]:= SolveAlways[a x + b == 0, x]
Out[1]= {{b -> 0, a -> 0}}

In[2]:= SolveAlways[(a + b) x + (a - b) y == 0, {x, y}]
Out[2]= {{a -> 0, b -> 0}}

In[3]:= SolveAlways[{a x + b == 0, c x + d == 0}, x]
Out[3]= {{b -> 0, a -> 0, d -> 0, c -> 0}}

In[4]:= SolveAlways[(a - b) x == 0, x]
Out[4]= {{b -> a}}
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
