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
