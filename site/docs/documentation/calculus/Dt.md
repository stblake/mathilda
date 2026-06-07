# Dt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Dt[f] gives the total derivative of f.
Dt[f, x] gives the total derivative of f with respect to x.
Dt[f, {x, n}] gives the nth total derivative.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Dt[y^2 + Sin[x]]
Out[1]= Cos[x] Dt[x] + 2 y Dt[y]

In[2]:= Dt[Pi + 3 + x y]
Out[2]= Dt[x] y + x Dt[y]

In[3]:= Dt[y^2, x]
Out[3]= 0

In[4]:= Dt[x^2, {x, 2}]
Out[4]= 2
```

## Implementation notes

- `Protected`, `ReadProtected`.
- Shares the elementary-function derivative table with `D`; the

**Attributes:** `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 2.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Dt[x y]
Out[1]= Dt[x] y + x Dt[y]
```

```mathematica
In[1]:= Dt[Sin[x]]
Out[1]= Cos[x] Dt[x]
```

```mathematica
In[1]:= Dt[Log[x]]
Out[1]= Dt[x]/x
```

```mathematica
In[1]:= Dt[a x, x]
Out[1]= a
```

### Notes

`Dt[f]` computes the total differential, treating every symbol as a potential independent variable and emitting `Dt[var]` factors for each one — so `Dt[x y]` gives the full product-rule expansion `Dt[x] y + x Dt[y]`. The two-argument form `Dt[f, x]` is the total derivative with respect to `x`, where other symbols are taken as constants unless they implicitly depend on `x`; `Dt[a x, x]` therefore returns `a`. Elementary functions differentiate through the chain rule with a residual `Dt[x]` factor, as in `Dt[Sin[x]]` and `Dt[Log[x]]`. Constants differentiate to `0` (`Dt[c, x]` gives `0`).
