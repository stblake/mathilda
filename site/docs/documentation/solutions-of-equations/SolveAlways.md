# SolveAlways

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SolveAlways[eqns, vars]`**

finds values of parameters appearing in eqns but not in vars such that eqns hold for every value of vars.

<details>
<summary>Notes</summary>

Equations may be Equal\[lhs, rhs\] or a List/And of such. Reduction: each lhs - rhs is treated as a polynomial in vars via CoefficientList; every coefficient must vanish, and the resulting system is passed to Solve with the remaining symbols (the parameters) as unknowns.  Returns {} when there are no parameters to solve for.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

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

### Applications (3)

```mathematica
In[5]:= SolveAlways[a x^2 + b x + c == 0, x]
Out[5]= {{c -> 0, b -> 0, a -> 0}}

In[6]:= SolveAlways[(a + b) x + (a - b - 2) == 0, x]
Out[6]= {{a -> 1, b -> -1}}

In[7]:= SolveAlways[p x^2 + q x + r == (x - 1)(x - 2), x]
Out[7]= {{r -> 2, q -> -3, p -> 1}}
```

## Algorithm

solvealways.c

Implementation of `SolveAlways[eqns, vars]`.

For each equation `lhs == rhs` we form the polynomial `p = lhs - rhs`, treat `p` as a polynomial in `vars` via `CoefficientList[p, vars]`,

```text
and require every coefficient to vanish.  The remaining symbols (those
```

appearing in `eqns` but not in `vars`) become the "parameters"; the collected coefficient equations are then passed to `Solve` with the parameters as unknowns.

Scope (v1):

```text
  - Equations may be a single `Equal[lhs, rhs]`, a `List[Equal, ...]`,
    or an `And[Equal, ...]`.
  - Variables may be a single symbol or a `List` of symbols.
  - Inequations (`Unequal`), disjunctions (`Or`), radicals, and
    `Series` strip are NOT handled here; they are deferred.
```

## Implementation notes

**Algorithm.** `builtin_solvealways` solves `SolveAlways[eqns, vars]` — find the parameter values making `eqns` hold for *all* values of `vars`. For each equation `lhs == rhs` it forms `p = lhs - rhs`, treats `p` as a polynomial in `vars` via `CoefficientList[p, vars]`, and requires every coefficient to vanish. The remaining symbols (appearing in `eqns` but not in `vars`) are the parameters; the collected coefficient equations are then handed to `Solve` with the parameters as the unknowns. Equations may be a single `Equal`, a `List` of `Equal`s, or an `And` of `Equal`s; variables a single symbol or a list of symbols.

**Data structures.** `Expr*`; coefficient extraction via `CoefficientList`, downstream solving delegated to the `Solve` builtin. Diagnostics use a one-shot hash-dedup pattern like `solve.c`.

**Limits.** v1 scope: inequations (`Unequal`), disjunctions (`Or`), radicals, and `Series` stripping are not handled and are deferred.

**Attributes:** `Protected`.

## References

**See also:** [Solve](../../solutions-of-equations/Solve/), [List](../../other-advanced/List/), [Unequal](../../comparisons/Unequal/), [Series](../../power-series/Series/), [CoefficientList](../../algebra/CoefficientList/), [Equal](../../comparisons/Equal/)

- Source: [`src/solvealways.c`](https://github.com/stblake/mathilda/blob/main/src/solvealways.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_solvealways.c`](https://github.com/stblake/mathilda/blob/main/tests/test_solvealways.c)

## Notes & additional examples

### Notes

`SolveAlways[eqns, vars]` finds the parameters (symbols in `eqns` but not in `vars`) for which the equations hold identically in `vars`. Each `lhs - rhs` is expanded as a polynomial in `vars` via `CoefficientList`; every coefficient is required to vanish, and the resulting system is handed to `Solve` for the parameters. It returns `{}` when there are no parameters to solve for.
