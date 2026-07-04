# Catalan

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Catalan
    is Catalan's constant, with numerical value ~= 0.915966.
Catalan is the sum over k >= 0 of (-1)^k (2 k + 1)^-2. It is a
mathematical constant: it has attributes Constant and Protected,
NumericQ[Catalan] is True, and D[Catalan, x] is 0. N[Catalan, prec]
evaluates it to any precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[Catalan] = {Constant,

**Attributes:** `Constant`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= N[Catalan]
Out[1]= 0.915966
```

```mathematica
In[1]:= N[Catalan, 40]
Out[1]= 0.91596559417721901505460351493238411077416
```

```mathematica
In[1]:= D[Catalan, x]
Out[1]= 0
```

```mathematica
In[1]:= NumericQ[Catalan]
Out[1]= True
```

```mathematica
In[1]:= N[8 Catalan, 30]
Out[1]= 7.32772475341775212043682811946
```

### Notes

`Catalan` is Catalan's constant G, the alternating sum over `k >= 0` of
`(-1)^k (2 k + 1)^-2`. It is a first-class symbolic constant: it carries the
`Constant` attribute (so `D[Catalan, x]` is `0`), satisfies `NumericQ`, and
`N[Catalan, prec]` evaluates it to any requested precision via an internal
fast-converging series. The 40-digit value matches the standard reference
value, and arithmetic combinations such as `8 Catalan` are held symbolically
until numericalized.
