# SolveAlways

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SolveAlways[eqns, vars]
    finds values of parameters appearing in eqns but not in vars
    such that eqns hold for every value of vars.
Equations may be Equal[lhs, rhs] or a List/And of such.
Reduction: each lhs - rhs is treated as a polynomial in vars
    via CoefficientList; every coefficient must vanish, and the
    resulting system is passed to Solve with the remaining
    symbols (the parameters) as unknowns.  Returns {} when there
    are no parameters to solve for.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SolveAlways[a x + b == 0, x]
Out[1]= {{b -> 0, a -> 0}}

In[2]:= SolveAlways[(a + b) x + (a - b) y == 0, {x, y}]
Out[2]= {{a -> 0, b -> 0}}

In[3]:= SolveAlways[{a x + b == 0, c x + d == 0}, x]
Out[3]= {{b -> 0, a -> 0, d -> 0, c -> 0}}

In[4]:= SolveAlways[(a - b) x == 0, x]
Out[4]= {{b -> a}}
```

## Implementation notes

**Algorithm.** `builtin_solvealways` solves `SolveAlways[eqns, vars]` — find the parameter values making `eqns` hold for *all* values of `vars`. For each equation `lhs == rhs` it forms `p = lhs - rhs`, treats `p` as a polynomial in `vars` via `CoefficientList[p, vars]`, and requires every coefficient to vanish. The remaining symbols (appearing in `eqns` but not in `vars`) are the parameters; the collected coefficient equations are then handed to `Solve` with the parameters as the unknowns. Equations may be a single `Equal`, a `List` of `Equal`s, or an `And` of `Equal`s; variables a single symbol or a list of symbols.

**Data structures.** `Expr*`; coefficient extraction via `CoefficientList`, downstream solving delegated to the `Solve` builtin. Diagnostics use a one-shot hash-dedup pattern like `solve.c`.

**Limits.** v1 scope: inequations (`Unequal`), disjunctions (`Or`), radicals, and `Series` stripping are not handled and are deferred.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/solvealways.c`](https://github.com/stblake/mathilda/blob/main/src/solvealways.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
