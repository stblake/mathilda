# Distribute

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Distribute[f[x1, x2, ...]]
    distributes f over Plus appearing in any of the xi, building the sum
    of f applied to every Cartesian-product selection of summands.
Distribute[expr, g]
    distributes over the head g instead of Plus.
Distribute[expr, g, f]
    performs the distribution only if the head of expr is f.
Distribute[expr, g, f, gp, fp]
    gives gp and fp in place of g and f respectively in the result.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Distribute[(a + b) . (x + y + z)]
Out[1]= Dot[a, x] + Dot[b, x] + Dot[a, y] + Dot[b, y] + Dot[a, z] + Dot[b, z]

In[2]:= Distribute[f[a + b, c + d + e]]
Out[2]= f[a, c] + f[b, c] + f[a, d] + f[b, d] + f[a, e] + f[b, e]

In[3]:= Distribute[(a + b + c) (u + v), Plus, Times]
Out[3]= a u + b u + c u + a v + b v + c v

In[4]:= Distribute[{{a, b}, {x, y, z}, {s, t}}, List]
Out[4]= {{a, x, s}, {a, x, t}, {a, y, s}, {a, y, t}, {a, z, s}, {a, z, t}, {b, x, s}, {b, x, t}, {b, y, s}, {b, y, t}, {b, z, s}, {b, z, t}}

In[5]:= Distribute[{{}, {a}}, {{}, {b}}, {{}, {c}}, List, List, List, Join]
Out[5]= Distribute[{{}, {a}}, {{}, {b}}, {{}, {c}}, List, List, List, Join]
```

## Implementation notes

- `Protected`.
- `Distribute` explicitly constructs the complete result of a distribution; `Expand`, on the other hand, builds up results iteratively, simplifying at each stage.
- For pure products, `Distribute` gives the same results as `Expand`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
