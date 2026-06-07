---
source: src/funcprog.c
---
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
