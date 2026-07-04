# FactorTermsList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FactorTermsList[poly]
    gives a list in which the first element is the overall numerical factor in poly, and the second element is the polynomial with the overall factor removed.
FactorTermsList[poly, {x1, x2, ...}]
    gives a list of factors of poly. The first element in the list is the overall numerical factor. The second element is a factor that does not depend on any of the xi. Subsequent elements are factors which depend on progressively more of the xi.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FactorTermsList[3 + 6 x + 3 x^2]
Out[1]= {3, 1 + 2 x + x^2}

In[2]:= FactorTermsList[14 x + 21 y + 35 x y + 63]
Out[2]= {7, 9 + 2 x + 3 y + 5 x y}

In[3]:= FactorTermsList[3 + 3 a + 6 a x + 6 x + 12 a x^2 + 12 x^2, x]
Out[3]= {3, 1 + a, 1 + 2 x + 4 x^2}

In[4]:= Times @@ FactorTermsList[14 x + 21 y + 35 x y + 63]
Out[4]= 7 (9 + 2 x + 3 y + 5 x y)
```

## Implementation notes

**Algorithm.** `builtin_factortermslist` (in `src/poly/facpoly_factorterms.inc`) is the list-returning sibling of `FactorTerms`. It calls the same engine `ft_compute_list` but returns the factor `List` directly instead of multiplying it out: `{c_0, c_1, …, c_k, residue}`, where `c_0` is the numerical content, `c_1..c_k` are the contents extracted with respect to progressively smaller subsets of the requested variables, and `residue` is the remaining polynomial part. Unlike `FactorTerms` it does not auto-thread (its result shape is already a list). See `FactorTerms` for the content-extraction algorithm.

**Data structures.** `Expr*` factors gathered into a single `List`; content computation uses the recursive multivariate-GCD helper `ft_content_wrt_set`.

- `Protected`.
- The product of the returned list always reproduces the input (after

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/poly/facpoly_factorterms.inc`](https://github.com/stblake/mathilda/blob/main/src/poly/facpoly_factorterms.inc)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= FactorTermsList[6 x^2 + 4 x]
Out[1]= {2, 2 x + 3 x^2}

In[2]:= FactorTermsList[2 x^2 + 4 x + 2]
Out[2]= {2, 1 + 2 x + x^2}
```

```mathematica
In[1]:= FactorTermsList[a b x^2 + a b c x, {x}]
Out[1]= {1, a b, c x + x^2}
```

```mathematica
In[1]:= FactorTermsList[12 x^3 y + 18 x^2 y, {x, y}]
Out[1]= {6, 1, y, 3 x^2 + 2 x^3}
```

### Notes

`FactorTermsList[poly]` returns `{overall numerical factor, polynomial with that factor removed}`. The numerical content (here 2) is pulled out, leaving the primitive part; it does not factor the remaining polynomial further. Given a variable list `{x1, ..., xn}`, it stratifies the polynomial: the first element is the numeric content, the second is the factor depending on none of the variables, and the later elements depend on progressively more of them. In the last example `12 x^3 y + 18 x^2 y` separates as `6 · 1 · y · (3 x^2 + 2 x^3)`, isolating the variable-free numeric content `6`, the `y`-only factor, and finally the `x`-dependent remainder.
