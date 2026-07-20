# PolynomialGCD

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PolynomialGCD[poly1, poly2, ...] gives the greatest common divisor of the polynomials.
Option Extension -> alpha computes the GCD over Q(alpha), where alpha is an
algebraic number recognised by qa_resolve_extension (Sqrt[c], c^(1/n), or I).
Default Extension -> None and Extension -> Automatic compute over the rationals,
treating any algebraic numbers in the input as independent variables.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PolynomialGCD[(1+x)^2(2+x)(4+x), (1+x)(2+x)(3+x)]
Out[1]= 2 + 3 x + x^2

In[2]:= PolynomialGCD[x^2+4x+4, x^2+2x+1]
Out[2]= 1

In[3]:= PolynomialGCD[x^2-1, x^3-1, x^4-1, x^5-1, x^6-1, x^7-1]
Out[3]= -1 + x

In[4]:= PolynomialGCD[x^2 - 2, x - Sqrt[2], Extension -> Sqrt[2]]
Out[4]= -Sqrt[2] + x

In[5]:= PolynomialGCD[x^3 - 2, x - 2^(1/3), Extension -> 2^(1/3)]
Out[5]= -2^(1/3) + x
```

## Implementation notes

**Algorithm.** `builtin_polynomialgcd` first strips an optional `Extension -> α` (or
`Extension -> Automatic`) and, when present, routes through the algebraic-number machinery
(`polynomialgcd_with_extension`, which lifts each input to a `QAUPoly` over `Q(α)` and folds
them with `qaupoly_gcd`; `qa_polynomialgcd_with_tower*` for multi-generator towers). Inexact
(floating) coefficients are force-rationalised, run through the exact algorithm, then
numericalised (`internal_rationalize_then_numericalize`).

The core path pre-processes each input with `decompose_to_bp` into a base/power list, peeling
off the integer content (numeric GCD of literal coefficients, including the integer content of
`Plus`-headed factors so it isn't double-counted) and any non-numeric factors common to every
argument. The remaining symbolic GCD is computed by `poly_gcd_internal`, which implements the
**recursive multivariate subresultant PRS**: it treats the last variable as main, splits each
operand into content (GCD of its coefficients, computed recursively in one fewer variable) and
primitive part, then reduces the primitive parts with `pseudo_rem` (a pseudo-remainder that
stays inside the coefficient ring, avoiding rationals) until the remainder is zero. The base
case (zero variables) is integer GCD via `my_number_gcd`. The result is `content_GCD ×
primitive_GCD`, normalised to a positive leading coefficient and expanded. Multi-argument GCD
folds left-to-right. A size budget (`max(input_size, 2000)` leaves) and a 50-iteration cap
guard against coefficient explosion over multi-radical rings; on overflow it conservatively
returns just the content GCD (always a valid divisor).

**Data structures.** Inputs are `Expr` trees; `BPList` holds the base/power decomposition;
`QAUPoly`/`QAExt`/`QATower` carry the algebraic-extension representation. Coefficients are
ordinary `Expr` subtrees, so coefficient GCDs recurse through the same machinery.

- `Protected`, `Listable`.
- Handles univariate and multivariate polynomials.
- Treats algebraic numbers (like `I`) as independent variables or constants seamlessly during complex arithmetic evaluations.
- Pre-extracts common factors before falling back to a full primitive Euclidean algorithm computation.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan): computes the GCD over `Q(alpha)` for `alpha` ∈ {`Sqrt[c]`, `c^(1/n)`, `I`} via lifting both inputs into the QAUPoly substrate (`src/poly/qaupoly.h`) and folding `qaupoly_gcd`. Extension support requires univariate inputs (after stripping the alpha-render symbol). Defaults `Extension -> None` and `Extension -> Automatic` work over the rationals and treat algebraic numbers as opaque variables. `Extension -> {alpha_1, ..., alpha_n}` (tower form) currently falls back to the no-extension path; tower-aware GCD is a Phase 0.5 follow-up.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- von zur Gathen & Gerhard, "Modern Computer Algebra" (3rd ed.), Ch. 6 & 11 (Euclidean and modular GCD).
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 7 (polynomial GCD computation).
- W. S. Brown, "On Euclid's Algorithm and the Computation of Polynomial Greatest Common Divisors", JACM 18(4), 1971.
- G. E. Collins, "Subresultants and Reduced Polynomial Remainder Sequences", JACM 14(1), 1967.
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= PolynomialGCD[x^2 - 1, x^2 + 2 x + 1]
Out[1]= 1 + x
```

```mathematica
In[1]:= PolynomialGCD[x^4 - 1, x^2 - 1]
Out[1]= -1 + x^2
```

```mathematica
In[1]:= PolynomialGCD[x^2 - 1, x - 1]
Out[1]= -1 + x
```

```mathematica
In[1]:= PolynomialGCD[x^3 - x, x^2 - x]
Out[1]= -x + x^2
```

```mathematica
In[1]:= PolynomialGCD[x^6 - 1, x^4 - 1, x^9 - 1]
Out[1]= -1 + x
```

```mathematica
In[1]:= PolynomialGCD[x^4 - 2, x^2 - Sqrt[2], Extension -> Sqrt[2]]
Out[1]= -Sqrt[2] + x^2
```

### Notes

`PolynomialGCD` returns the greatest common divisor of its polynomial
arguments, here over the rationals. The result is the highest-degree common
factor — for `x^2 - 1` and `(x+1)^2` that shared factor is `1 + x`, and for
`x^4 - 1` and `x^2 - 1` it is the full `x^2 - 1`. The output is normalized in
canonical term order and is not forced monic, so a shared `x` factor surfaces
as `-x + x^2`. The `Extension -> alpha` option computes the GCD over `Q(alpha)`
for `Sqrt[c]`, `c^(1/n)`, or `I`; the default treats any algebraic numbers as
independent variables.
