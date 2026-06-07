# D

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
D[f, x] gives the partial derivative of f with respect to x.
D[f, {x, n}] gives the nth partial derivative.
D[f, x, y, ...] gives the mixed derivative.
D[f, x, NonConstants -> {y, ...}] treats the listed symbols as implicit functions of x.
Distributes over Equal: D[a == b, x] gives D[a, x] == D[b, x].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= D[x^3, x]
Out[1]= 3 x^2

In[2]:= D[Sin[x^2], x]
Out[2]= 2 x Cos[x^2]

In[3]:= D[Sin[a x], {x, 3}]
Out[3]= -a^3 Cos[a x]

In[4]:= D[f[g[x]], x]
Out[4]= Derivative[1][g][x] Derivative[1][f][g[x]]

In[5]:= D[f[x, y], x]
Out[5]= Derivative[1, 0][f][x, y]

In[6]:= D[Derivative[2][f][x], x]
Out[6]= Derivative[3][f][x]

In[7]:= D[{x, x^2, Sin[x]}, x]
Out[7]= {1, 2 x, Cos[x]}

In[8]:= D[Log[b, x], x]
Out[8]= 1/(Log[b] x)
```

## Implementation notes

- `Protected`, `ReadProtected`.
- Recognises the elementary heads `Plus`, `Times`, `Power`, `Sqrt`,

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
In[1]:= D[x^n, x]
Out[1]= n x^(-1 + n)
```

```mathematica
In[1]:= D[Exp[a x], x]
Out[1]= a E^(a x)
```

```mathematica
In[1]:= D[x^2 y, {x, 2}]
Out[1]= 2 y
```

```mathematica
In[1]:= D[f[g[x]], x]
Out[1]= Derivative[1][g][x] Derivative[1][f][g[x]]
```

### Notes

`D` is driven by the pattern-based derivative rules in `src/internal/deriv.m`, implementing the chain, product, and quotient rules together with the elementary-function table. Differentiating an unknown function head produces a `Derivative[n][f]` operator rather than evaluating further, so the chain rule on `f[g[x]]` returns a product of such operators. The `{x, n}` form takes the `n`th derivative and treats other symbols as constants by default; use `NonConstants -> {...}` to mark implicit dependencies.
