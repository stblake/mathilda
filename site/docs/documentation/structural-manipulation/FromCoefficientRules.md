# FromCoefficientRules

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FromCoefficientRules[{expvec -> coeff, ...}, {x1, x2, ...}] reconstructs the polynomial`**

**`Sum[coeff * x1^e1 * x2^e2 * ..., over rules]. The exponent vectors must match the number of`**

<details>
<summary>Notes</summary>

variables. Inverse of CoefficientRules.

</details>

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= FromCoefficientRules[{{2, 0} -> a, {1, 1} -> b, {0, 2} -> c}, {x, y}]
Out[1]= a x^2 + b x y + c y^2
```

## Algorithm

monomials.c — MonomialList, CoefficientRules, FromCoefficientRules

A sparse {exponent-vector -> coefficient} view of a multivariate polynomial, plus the list of its monomials and the inverse reconstruction. All three heads share one core:

```text
  1. resolve the variable list (explicit, `All`, or default Variables[poly]);
  2. reduce coefficients modulo an optional `Modulus -> m`;
  3. Expand and split into additive terms;
  4. decompose each term into (integer exponent vector, coefficient) w.r.t.
     the variables;
  5. merge like monomials and drop zero coefficients;
  6. sort by a monomial order (six named orders + explicit weight matrix).
```

MonomialList and CoefficientRules differ only in how each (expvec, coeff) term is rendered. FromCoefficientRules is the inverse of CoefficientRules.

Every named order is a special case of "descending lexicographic order of the weighted exponent vectors w.v" (Wolfram's own model), so a single weight matrix drives one comparator. For k variables (e_i = i-th unit row, deg = all-ones row), greatest monomial first:

```text
  Lexicographic                    e_0, e_1, ..., e_{k-1}          (default)
  NegativeLexicographic           -e_0, ..., -e_{k-1}   (= Sort, ascending)
  DegreeLexicographic              deg, e_0, ..., e_{k-2}
  DegreeReverseLexicographic       deg, -e_{k-1}, ..., -e_1
  NegativeDegreeLexicographic     -deg, e_0, ..., e_{k-2}
  NegativeDegreeReverseLexicographic  -deg, -e_{k-1}, ..., -e_1
```

These reproduce Wolfram's explicit-matrix spellings exactly (e.g. DegreeLexicographic on {x,y} is {{1,1},{1,0}}; DegreeReverseLexicographic on {x,y,z} is {{1,1,1},{0,0,-1},{0,-1,0}}).

## Implementation notes

- `Protected`.
- Builds `Sum[coeff * var1^e1 * var2^e2 * ..., over rules]`; the exponent vectors must
  match the number of variables.
- Inverse of `CoefficientRules`: `FromCoefficientRules[CoefficientRules[poly, vars], vars]`
  is `Expand[poly]`.

**Attributes:** `Protected`.

## References

**See also:** [CoefficientRules](../../structural-manipulation/CoefficientRules/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_coefficient_rules.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coefficient_rules.c)
