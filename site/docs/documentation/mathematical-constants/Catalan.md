# Catalan

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Catalan`**

is Catalan's constant, with numerical value ~= 0.915966.

**`NumericQ[Catalan] is True, and D[Catalan, x] is 0. N[Catalan, prec]`**

<details>
<summary>Notes</summary>

Catalan is the sum over k \>= 0 of (-1)^k (2 k + 1)^-2. It is a mathematical constant: it has attributes Constant and Protected, evaluates it to any precision.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= N[Catalan]
Out[1]= 0.915966

In[2]:= N[Catalan, 40]
Out[2]= 0.91596559417721901505460351493238411077416

In[3]:= D[Catalan, x]
Out[3]= 0

In[4]:= NumericQ[Catalan]
Out[4]= True

In[5]:= N[8 Catalan, 30]
Out[5]= 7.32772475341775212043682811946
```

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[Catalan] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[Catalan]` is `True` and
  `D[Catalan, x] = 0`.
- `N[Catalan]` gives the machine value `0.915966`; `N[Catalan, prec]` gives any
  precision (MPFR `mpfr_const_catalan`), e.g.
  `N[Catalan, 50] = 0.915965594177219015054603514932384110774149374281673`.

**Attributes:** `Constant`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Notes

`Catalan` is Catalan's constant G, the alternating sum over `k >= 0` of
`(-1)^k (2 k + 1)^-2`. It is a first-class symbolic constant: it carries the
`Constant` attribute (so `D[Catalan, x]` is `0`), satisfies `NumericQ`, and
`N[Catalan, prec]` evaluates it to any requested precision via an internal
fast-converging series. The 40-digit value matches the standard reference
value, and arithmetic combinations such as `8 Catalan` are held symbolically
until numericalized.
