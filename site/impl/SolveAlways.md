---
source: src/solvealways.c
---
**Algorithm.** `builtin_solvealways` solves `SolveAlways[eqns, vars]` — find the parameter values making `eqns` hold for *all* values of `vars`. For each equation `lhs == rhs` it forms `p = lhs - rhs`, treats `p` as a polynomial in `vars` via `CoefficientList[p, vars]`, and requires every coefficient to vanish. The remaining symbols (appearing in `eqns` but not in `vars`) are the parameters; the collected coefficient equations are then handed to `Solve` with the parameters as the unknowns. Equations may be a single `Equal`, a `List` of `Equal`s, or an `And` of `Equal`s; variables a single symbol or a list of symbols.

**Data structures.** `Expr*`; coefficient extraction via `CoefficientList`, downstream solving delegated to the `Solve` builtin. Diagnostics use a one-shot hash-dedup pattern like `solve.c`.

**Limits.** v1 scope: inequations (`Unequal`), disjunctions (`Or`), radicals, and `Series` stripping are not handled and are deferred.
