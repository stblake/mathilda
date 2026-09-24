# MonomialList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MonomialList[poly] gives the list of monomials of poly, using Variables[poly].`**

**`MonomialList[poly, {x1, x2, ...}] uses the given variables; MonomialList[poly, vars, order]`**

<details>
<summary>Notes</summary>

sorts by order. order is "Lexicographic" (default), "DegreeLexicographic", "DegreeReverseLexicographic", "NegativeLexicographic", "NegativeDegreeLexicographic", "NegativeDegreeReverseLexicographic", or an explicit weight matrix. Modulus -\> m reduces coefficients modulo m. vars may be All (equivalent to Variables\[poly\]).

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= MonomialList[(x + y)^3]
Out[1]= {x^3, 3 x^2 y, 3 x y^2, y^3}

In[2]:= MonomialList[x^2 y^2 + x^3, {x, y}, "DegreeLexicographic"]
Out[2]= {x^2 y^2, x^3}
```

### Options (1)

```mathematica
In[3]:= MonomialList[(x + 1)^5, x, Modulus -> 2]
Out[3]= {x^5, x^4, x, 1}
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
- `MonomialList[poly]` is equivalent to `MonomialList[poly, Variables[poly]]`.
- Works whether or not `poly` is given in expanded form (it is expanded internally).
- `order` is `"Lexicographic"` (default), `"DegreeLexicographic"`,
  `"DegreeReverseLexicographic"`, `"NegativeLexicographic"`,
  `"NegativeDegreeLexicographic"`, `"NegativeDegreeReverseLexicographic"`, or an
  explicit weight matrix `w` (monomials ranked by the lexicographic order of `w.v`).
- `"NegativeLexicographic"` is `Sort` of the exponent vectors; `"Lexicographic"` is
  its reverse.
- `Modulus -> m` reduces coefficients modulo `m` (dropping monomials whose
  coefficient vanishes).
- `vars` may be `All` (equivalent to `Variables[poly]`).
- `Plus @@ MonomialList[poly, vars]` reconstructs the expanded polynomial.

**Attributes:** `Protected`.

## References

**See also:** [Sort](../../data-structures/Sort/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_coefficient_rules.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coefficient_rules.c)
