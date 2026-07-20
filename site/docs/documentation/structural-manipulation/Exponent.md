# Exponent

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Exponent[expr, form] gives the maximum power with which form appears in the
expanded form of expr.
Exponent[expr, form, h] applies h to the set of exponents (default Max).
form may be a symbol, a kernel, or a product of terms; expr need not be
expanded.  Exponent is purely syntactic (no zero-coefficient recognition).
Exponent[0, x] is -Infinity.  Exponent[expr, {f1, f2, ...}] gives the list of
exponents for each fi.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Exponent[1 + x^2 + a x^3, x]
Out[1]= 3

In[2]:= Exponent[0, x]
Out[2]= -Infinity

In[3]:= Exponent[x^(n+1) + 2 Sqrt[x] + 1, x]
Out[3]= Max[1/2, 1 + n]

In[4]:= Exponent[(x^2+1)^3 - 1, x, Min]
Out[4]= 2

In[5]:= Exponent[1 + x^2 + a x^3, x, List]
Out[5]= {0, 2, 3}
```

## Implementation notes

- `Listable`, `Protected`.
- The default aggregator is `h = Max`. `Exponent[expr, form, Min]` gives the lowest power; `Exponent[expr, form, List]` gives the sorted, de-duplicated set of exponents.
- `form` may be a symbol, a kernel (e.g. `Sin[x]`), or a product of terms.
- Works whether or not `expr` is explicitly given in expanded form (it expands internally).
- Purely syntactic: it does not attempt to recognise zero coefficients.
- Exponents may be rational numbers or symbolic expressions.
- `Exponent[0, x]` is `-Infinity` (empty exponent set, `h = Max`).
- The `Listable` attribute makes `Exponent[expr, {form1, form2, ...}]` give the list of exponents for each `formi`.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
