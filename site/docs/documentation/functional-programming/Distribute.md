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

**Algorithm.** `builtin_distribute` performs the full combinatorial expansion of
a head `f` over a "sum-like" head `g` (default `Plus`), e.g. `Distribute[f[a+b,
c+d]]` → `f[a,c]+f[a,d]+f[b,c]+f[b,d]`. It accepts up to five arguments `expr, g,
f, gp, fp`, where `gp`/`fp` are the heads used to *rebuild* the result (default
`gp=g`, `fp=f`). It first checks `Head[expr] == f`; if not, `expr` is returned
unchanged.

Each argument of `expr` becomes a *component list*: arguments whose head equals
`g` contribute all their summands, every other argument contributes itself as a
single-element component. If no argument is a `g`-expression there is nothing to
distribute and `expr` is returned. Otherwise `distribute_recursive` forms the
Cartesian product of the component lists, emitting one `fp[tuple...]` term per
combination; these terms are collected (with a doubling-capacity buffer) and
wrapped in `gp[...]`. The result is run through `evaluate()` because `gp`/`fp`
may be builtins (`Plus`, `Times`, ...).

**Complexity.** Output size is the product of the component-list lengths, so the
term count grows multiplicatively in the number of summands per slot — the
expansion is genuinely combinatorial.

**Data structures.** Component lists are arrays of borrowed `Expr*` (sub-args of
`expr`); the recursion threads a `current_tuple` scratch array and a growable
`results` array, both freed after the wrapper `gp[...]` is built.

- `Protected`.
- `Distribute` explicitly constructs the complete result of a distribution; `Expand`, on the other hand, builds up results iteratively, simplifying at each stage.
- For pure products, `Distribute` gives the same results as `Expand`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Distribute[(a + b)(c + d)]
Out[1]= a c + b c + a d + b d
```

```mathematica
In[1]:= Distribute[f[a + b, c + d]]
Out[1]= f[a, c] + f[b, c] + f[a, d] + f[b, d]
```

```mathematica
In[1]:= Distribute[And[a, Or[b, c]], Or, And]
Out[1]= a && b || a && c
```

```mathematica
In[1]:= Distribute[f[a + b + c], Plus, f, Plus, g]
Out[1]= g[a] + g[b] + g[c]
```

### Notes

`Distribute[f[x1, x2, ...]]` forms the sum of `f` applied to every term in the
Cartesian product of the summands found in the `xi`. On `Times` it reproduces
ordinary expansion (`(a+b)(c+d)`), but the mechanism is general: `f[a+b, c+d]`
distributes an arbitrary head over its `Plus` arguments. The two extra-argument
form `Distribute[expr, g]` distributes over the head `g` instead of `Plus`; the
three-argument form `Distribute[expr, g, f]` only acts when the head of `expr` is
`f`. The third example distributes `And` over `Or` — the logical distributive law,
turning a clause into disjunctive form. The five-argument form
`Distribute[expr, g, f, gp, fp]` substitutes `gp` and `fp` for `g` and `f` in the
result, here rewriting the distributed terms under a new function head `g`.
