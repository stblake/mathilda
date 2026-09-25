# PolynomialReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PolynomialReduce[poly, {p1, ..., pn}, {x1, ..., xk}] gives {{a1, ..., an}, b} with`**

<details>
<summary>Notes</summary>

a1 p1 + ... + an pn + b == poly and b minimal: no term of b is divisible by any leading term of the pi. PolynomialReduce\[poly, polys\] uses Variables. Options as for GroebnerBasis: MonomialOrder (Lexicographic default), CoefficientDomain (RationalFunctions default; Rationals), Modulus -\> p (over GF(p)), ParameterVariables. Free symbols outside the variable list are coefficient-field parameters. If the pi are a Groebner basis, b is the unique normal form.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= PolynomialReduce[x^3 + y^3, {x^2 - y^2 - 1, x + 2 y - 7}, {x, y}]
Out[1]= {{x, 1 + y^2}, 7 - 2 y + 7 y^2 - y^3}

In[2]:= f == First[%].{x^2 - y^2 - 1, x + 2 y - 7} + Last[%] // Expand
Out[2]= f == -1 + Dot[-1, {-1 + x^2 - y^2, -7 + x + 2 y}]

In[3]:= PolynomialReduce[2 x^3 + y^3 + 3 y, {x^2 + y^2 - 1, x y - 2}, {x, y}]
Out[3]= {{2 x, -2 y}, 2 x - y + y^3}

In[4]:= gb = GroebnerBasis[{x^2 - y^3 - 5, y^2 - x^3 + 7}, {x, y}]; PolynomialReduce[x^6 - 14 x^3 - y^4 + 49, gb, {x, y}][[2]]
Out[4]= 0

In[5]:= PolynomialReduce[x^5 + (x + y)^2, {x^2 - y^3 - 5, y^3 - x^3 + 7}, {y, x}]
Out[5]= {{0, 0}, x^2 + x^5 + 2 x y + y^2}
```

### Options (2)

```mathematica
In[6]:= gb1 = GroebnerBasis[{a x^2 + 5 y - 1, 2 x + x y - y^2}, {x, y}, CoefficientDomain -> RationalFunctions]; PolynomialReduce[a^2 x - x y + y^2 - 3, gb1, {x, y}][[1]]
Out[6]= {1/4/a, 1/4 a - 1/4 y/a}

In[7]:= PolynomialReduce[x^3, {x^2 + 1}, {x}, Modulus -> 7]
Out[7]= {{x}, 6 x}
```

## Options & behaviour

**Options** (as for `GroebnerBasis`):

| Option               | Default              | Supported                                          |
| -------------------- | -------------------- | -------------------------------------------------- |
| `MonomialOrder`      | `Lexicographic`      | `Lexicographic`, `DegreeLexicographic`, `DegreeReverseLexicographic`, or an integer weight matrix. Non-well-founded (`Negative*`) orders fall back to `Lexicographic`. |
| `CoefficientDomain`  | `RationalFunctions`  | `RationalFunctions` (the default; over the field `Q(params)`, so cofactors may carry parameter denominators), `Rationals` (over `Q`; requires every symbol to be a variable). |
| `Modulus`            | `0`                  | `Modulus -> p` (prime `p` in `[2, 2^31)`) reduces over `GF(p)` (`gfp_divmod`, `Lexicographic` / `DegreeReverseLexicographic`; not combined with parameters). |
| `ParameterVariables` | `{}`                 | Symbol or list of symbols treated as coefficient-field parameters. |
| `Tolerance`          | `0`                  | Only `0` is supported. |

**Deferred** (accepted, then the head declines with a
`PolynomialReduce::nimpl` note): `CoefficientDomain -> Integers` and
`InexactNumbers`, and a nonzero `Tolerance`.

## Algorithm

polynomialreduce.c

```text
`PolynomialReduce[poly, {p1, ..., pn}, {x1, ..., xk}]` -- multivariate
polynomial division with cofactors.  Returns {{a1, ..., an}, b} such that

    a1 p1 + a2 p2 + ... + an pn + b == poly   (exactly)
```

and b is fully reduced: no term of b is divisible by any leading term of the

```text
pi under the chosen monomial order.  This is the normal-form division that
```

sits at the heart of Buchberger's algorithm, exposed as a user builtin; it is the multivariate-division sibling of GroebnerBasis and shares its options.

Three coefficient domains are implemented:

```text
  - Rationals over Q (and RationalFunctions with no parameters, which is the
    same field): the exact GBPoly engine gb_divmod (groebner.c), with a FLINT
    fast path (fmpq_mpoly_divrem_ideal) for the Lexicographic case on a
    FLINT-enabled build.  gb_divmod is the authoritative, Wolfram-matching
    engine.
  - RationalFunctions over Q(params): when free symbols outside the variable
    list appear, they are coefficient-field elements, so the cofactors can
    carry denominators in the parameters (e.g. y/(4 a)).  The GBPoly engine
    (rational coefficients only) cannot express those; this case runs an
    Expr-coefficient division whose coefficient arithmetic is Together/Cancel
    over Q(params) -- exact, no numeric oracle.
  - Modulus -> p (GF(p)): via the gbmod.c gfp_divmod engine.
```

CoefficientDomain -> Integers and InexactNumbers (+ Tolerance) are distinct ring/numeric engines and are DEFERRED: the option is accepted but the head declines (returns NULL) with a note, so behaviour is never silently wrong.

```text
Attributes: Protected.  Options: MonomialOrder, CoefficientDomain, Modulus,
```

ParameterVariables, Tolerance (defaults registered in options_builtin.c).

PolynomialReduce is a symbolic/structural head (it returns lists of symbolic polynomials, not element-wise machine numbers over a buffer), so it is genuinely exempt from the packed/NDArray and Compile[] numeric surfaces -- there is nothing element-wise to lower.

## Implementation notes

- `Protected`.
- The engine is `gb_divmod` (`src/poly/groebner.c`), the same normal-form
  division Buchberger uses; over Q with `Lexicographic` order it dispatches
  to FLINT's `fmpq_mpoly_divrem_ideal` when FLINT is present.
- Shares `GroebnerBasis`'s monomial orders (including `DegreeLexicographic`
  via the shared named-order weight matrix) and its parameter discovery: a
  free symbol outside the variable list is a coefficient-field parameter.
- The quotients and the remainder depend on the monomial order and on the
  order of the `pi` (as in `GroebnerBasis` / the multivariate division
  algorithm).
- Inputs are `Expand`ed and (inexact leaves) rationalised before division.

**Attributes:** `Protected`.

## References

**See also:** [Variables](../../algebra/Variables/), [GroebnerBasis](../../algebra/GroebnerBasis/), [Expand](../../algebra/Expand/), [Modulus](../../other-advanced/Modulus/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_polynomialreduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_polynomialreduce.c)
- Tests: [`tests/test_risch_field.c`](https://github.com/stblake/mathilda/blob/main/tests/test_risch_field.c)
- Tests: [`tests/test_risch_hypertangent.c`](https://github.com/stblake/mathilda/blob/main/tests/test_risch_hypertangent.c)
