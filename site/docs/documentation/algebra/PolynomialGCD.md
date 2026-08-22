# PolynomialGCD

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PolynomialGCD[poly1, poly2, ...] gives the greatest common divisor of the polynomials.`**

<details>
<summary>Notes</summary>

Option Extension -\> alpha computes the GCD over Q(alpha), where alpha is an algebraic number recognised by qa\_resolve\_extension (Sqrt\[c\], c^(1/n), or I). Default Extension -\> None and Extension -\> Automatic compute over the rationals, treating any algebraic numbers in the input as independent variables.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= PolynomialGCD[(1+x)^2(2+x)(4+x), (1+x)(2+x)(3+x)]
Out[1]= 2 + 3 x + x^2

In[2]:= PolynomialGCD[x^2+4x+4, x^2+2x+1]
Out[2]= 1

In[3]:= PolynomialGCD[x^2-1, x^3-1, x^4-1, x^5-1, x^6-1, x^7-1]
Out[3]= -1 + x
```

### Options (2)

```mathematica
In[4]:= PolynomialGCD[x^2 - 2, x - Sqrt[2], Extension -> Sqrt[2]]
Out[4]= -Sqrt[2] + x

In[5]:= PolynomialGCD[x^3 - 2, x - 2^(1/3), Extension -> 2^(1/3)]
Out[5]= -2^(1/3) + x
```

### Applications (6)

```mathematica
In[6]:= PolynomialGCD[x^2 - 1, x^2 + 2 x + 1]
Out[6]= 1 + x

In[7]:= PolynomialGCD[x^4 - 1, x^2 - 1]
Out[7]= -1 + x^2

In[8]:= PolynomialGCD[x^2 - 1, x - 1]
Out[8]= -1 + x

In[9]:= PolynomialGCD[x^3 - x, x^2 - x]
Out[9]= -x + x^2

In[10]:= PolynomialGCD[x^6 - 1, x^4 - 1, x^9 - 1]
Out[10]= -1 + x

In[11]:= PolynomialGCD[x^4 - 2, x^2 - Sqrt[2], Extension -> Sqrt[2]]
Out[11]= -Sqrt[2] + x^2
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Discriminant of deg 20 | 2.51 s | 0.068 s | 0.182 s |
| Expand (1+x)^400 | 0.434 s | 0.107 s | 0.003 s |
| Cancel deg-60 over deg-58 | 0.337 s | 0.569 s | 7.37 s |
| PolynomialGCD, coprime deg 40 | 0.252 s | 0.087 s | 0.334 s |
| PolynomialGCD, shared deg-20 factor | 0.079 s | 0.063 s | 0.764 s |
| PolynomialQuotient deg 60 / deg 20 | 0.063 s | 0.209 s | 0.945 s |

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

## References

**See also:** [I](../../mathematical-constants/I/)

- von zur Gathen & Gerhard, "Modern Computer Algebra" (3rd ed.), Ch. 6 & 11 (Euclidean and modular GCD).
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 7 (polynomial GCD computation).
- W. S. Brown, "On Euclid's Algorithm and the Computation of Polynomial Greatest Common Divisors", JACM 18(4), 1971.
- G. E. Collins, "Subresultants and Reduced Polynomial Remainder Sequences", JACM 14(1), 1967.
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_eval_timestamps.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval_timestamps.c)
- Tests: [`tests/test_expr_sharing.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expr_sharing.c)
- Tests: [`tests/test_extension_auto_builtins.c`](https://github.com/stblake/mathilda/blob/main/tests/test_extension_auto_builtins.c)
- Tests: [`tests/test_extension_options.c`](https://github.com/stblake/mathilda/blob/main/tests/test_extension_options.c)

## Notes & additional examples

### Notes

`PolynomialGCD` returns the greatest common divisor of its polynomial
arguments, here over the rationals. The result is the highest-degree common
factor — for `x^2 - 1` and `(x+1)^2` that shared factor is `1 + x`, and for
`x^4 - 1` and `x^2 - 1` it is the full `x^2 - 1`. The output is normalized in
canonical term order and is not forced monic, so a shared `x` factor surfaces
as `-x + x^2`. The `Extension -> alpha` option computes the GCD over `Q(alpha)`
for `Sqrt[c]`, `c^(1/n)`, or `I`; the default treats any algebraic numbers as
independent variables.
