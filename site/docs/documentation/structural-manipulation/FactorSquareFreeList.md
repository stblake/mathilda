# FactorSquareFreeList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FactorSquareFreeList[poly] gives a list of the square-free factors of poly`**

<details>
<summary>Notes</summary>

together with their multiplicities, as {factor, exponent} pairs.  A thin wrapper over FactorSquareFree: it splits that product form into pairs and forwards the Extension option verbatim.  The first element is always the overall numerical factor {c, 1} (it is {1, 1} when there is none).

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= FactorSquareFreeList[x^5 - x^3 - x^2 + 1]
Out[1]= {{1, 1}, {-1 + x, 2}, {1 + 2 x + 2 x^2 + x^3, 1}}

In[2]:= FactorSquareFreeList[2 x^3 + 2 x^2 - 2 x - 2]
Out[2]= {{2, 1}, {-1 + x, 1}, {1 + x, 2}}
```

### Options (1)

```mathematica
In[3]:= FactorSquareFreeList[x^2 - 2 Sqrt[2] x + 2, Extension -> Sqrt[2]]
Out[3]= {{1, 1}, {-Sqrt[2] + x, 2}}
```

## Algorithm

facpoly_list.c -- FactorList[poly] / FactorList[poly, opts].

A thin wrapper over Factor: factor `poly` via `Factor[poly, opts...]` (all options are forwarded verbatim -- GaussianIntegers, Extension, ...), then split the resulting product into {factor, exponent} pairs.

```text
  FactorList[x^2 - 1]              -> {{1, 1}, {-1 + x, 1}, {1 + x, 1}}
  FactorList[2 x^3 + 2 x^2 - ...]  -> {{2, 1}, {-1 + x, 1}, {1 + x, 2}}
  FactorList[x^4 - 2, Extension -> Sqrt[2]]
                                   -> {{1, 1}, {Sqrt[2] + x^2, 1}, {-Sqrt[2] + x^2, 1}}
```

The first element is always the overall numerical factor {c, 1} (which is

```text
{1, 1} when there is no numerical factor).  Denominator factors of a
```

rational function carry negative exponents.

Parsing rules on the Factor output R (a Times, a bare factor, or a number):

```text
  - a number literal (Integer / Rational / Real / Complex / ...) multiplies
    into the overall numerical factor `c`;
  - Power[base, e] with an *integer* e is the pair {base, e} (this is a
    factor raised to a multiplicity, positive or negative);
  - anything else -- including Power[base, 1/2] = Sqrt[base], which is an
    irreducible factor in its own right -- is the pair {factor, 1}.
```

## Implementation notes

- `Listable`, `Protected`.
- A thin wrapper over `FactorSquareFree`: it decomposes via `FactorSquareFree[poly, opts...]` (the `Extension` option is forwarded verbatim) and splits the product into `{factor, exponent}` pairs.  The relationship to `FactorSquareFree` is exactly `FactorList`'s to `Factor`.
- The first element is always the overall numerical factor `{c, 1}` — it is `{1, 1}` when there is no numerical factor.
- Each square-free factor keeps its multiplicity; unlike `FactorList`, the factors are *not* further reduced to irreducibles (`x^2 - 1` stays whole rather than splitting into `(-1 + x)(1 + x)`).
- `Times @@ Power @@@ FactorSquareFreeList[poly]` reconstructs `poly` (up to `FactorSquareFree`'s normal form).

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [FactorSquareFree](../../algebra/FactorSquareFree/), [FactorList](../../structural-manipulation/FactorList/), [Factor](../../algebra/Factor/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_factorsquarefreelist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_factorsquarefreelist.c)
