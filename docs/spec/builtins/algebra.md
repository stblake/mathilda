# Algebra

Algebraic manipulation of rational and polynomial expressions: rational-expression normal forms (`Numerator`, `Denominator`, `Cancel`, `Together`, `Apart`), resultant-based tools, and minimal and irreducible-polynomial predicates. The equation solver (`Solve`, `Root`, `ToRadicals`, ...) lives in [`solutions-of-equations.md`](solutions-of-equations.md).

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

In[7]:= Together[1/(x - k^(1/3)) + 1/(x + k^(1/3))]
Out[7]= (2 x)/(-k^(2/3) + x^2)
```

- **Cube-and-higher-root towers** (2026-07-04): a `Plus` of fractions over
  radicals `Power[base, p/q]` of any index `q ≥ 2` whose radicand carries a free
  symbol (e.g. `k^(1/3)`) is combined over a common denominator via FLINT
  (`flint_algebraic_field_together`), with the radicals treated as free kernels —
  WL-faithful: radicals stay in the denominator, no rationalisation. This
  generalises the sqrt-only combiner to any index. The same path serves `Cancel`;
  single-fraction relation-dependent cancellations (`(x^3-k)/(x-k^(1/3))`) still
  go through the relation-aware GCD path.

## RootReduce

Canonicalises an algebraic expression. `RootReduce` dispatches between two
rigorous FLINT engines by the shape of the input:

- `RootReduce[expr]`
- `RootReduce[expr, Method -> "Automatic" | "Recursive" | "NumberField"]`

**Constant algebraic numbers** (no free symbol) — integers, rationals, radicals,
roots of unity, the imaginary unit and `Root[]` objects combined by `+,-,*,/,^`
— are canonicalised via FLINT `qqbar` (`src/poly/flint_qqbar.c`) to a single
representative: a **rational number**, a **quadratic radical expression**, or a
`Root[Function[minpoly&], k]` object (degree ≥ 3). The `Root` index `k` follows
the Wolfram ordering (real roots ascending first, then non-real roots).

**Algebraic functions** over a tower `Q(params)(radicals)` — radicals whose
radicand carries a free variable (e.g. the Goursat `k^(1/3)` towers) — have
their denominator **rationalised** by `flint_algebraic_field_canonical`
(`src/poly/flint_bridge.c`): the tower is a finite-dimensional `Q(params)`-vector
space, multiplication by the denominator is a linear map (products reduced modulo
the generators' minimal-polynomial ideal via `fmpz_mpoly_divrem_ideal`), and
inverting it with a solve over `Q(params)` (`gr_mat` over `fmpz_mpoly_q`)
rationalises the denominator to `Norm_{K/Q(params)}(D)`. No numeric zero oracle.

**Features**:
- `Protected`, `Listable`. Threads over lists, and over equations, inequalities
  and logic functions (`Equal`, `Unequal`, `Less`, `And`, ...); for
  (in)equalities of constant algebraic numbers it decides the relation exactly
  via `qqbar`.
- `Method`: `"Recursive"`/`"Automatic"` fold `qqbar` arithmetic bottom-up;
  `"NumberField"` re-expresses the value through a single primitive element of a
  common number field (`qqbar_express_in_field`). All three yield the identical
  canonical result. A `Root[]` object of degree ≤ 2 (or degree 1) auto-reduces
  to a quadratic radical / rational.
- One positional argument is required; other arg counts emit `RootReduce::argx`.
  An unknown `Method` emits `RootReduce::mtd`. Idempotent.

```mathematica
In[1]:= RootReduce[Sqrt[2] + Sqrt[3]]
Out[1]= Root[1 - 10 #1^2 + #1^4 &, 4]

In[2]:= RootReduce[(Sqrt[18] + Sqrt[27]) / Sqrt[5 + 2 Sqrt[6]]]
Out[2]= 3

In[3]:= RootReduce[1/(1 + Sqrt[2])]
Out[3]= -1 + Sqrt[2]

In[4]:= RootReduce[1/(1 + 2^(1/3) + 2^(2/3))]
Out[4]= Root[-1 + 3 #1 + 3 #1^2 + #1^3 &, 1]

In[5]:= RootReduce[Sqrt[2] + Sqrt[3] + Sqrt[5] ==
          Sqrt[10 + 2 Sqrt[15] + 4 Sqrt[4 + Sqrt[15]]]]
Out[5]= True

In[6]:= RootReduce[1/(1 + k^(1/3))]        (* parametric tower *)
Out[6]= (1 - k^(1/3) + k^(2/3)) / (1 + k)
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

