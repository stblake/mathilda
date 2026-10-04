# Activate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Activate[expr] reactivates every Inactive[h] in expr (replacing Inactive[h] by h) and re-evaluates, so an inactive integral Inactive[Integrate][g,x] becomes Integrate[g,x] and evaluates.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Activate[Inactive[Integrate][2 y, y]]
Out[1]= y^2
```

### Applications (2)

Reactivates the inert head, so the sum is taken

```mathematica
In[2]:= Activate[Inactive[Plus][2, 3]]
Out[2]= 5
```

An inert integral becomes a real one and evaluates

```mathematica
In[3]:= Activate[Inactive[Integrate][x, x]]
Out[3]= 1/2 x^2
```

## Implementation notes

**Algorithm.** `builtin_activate` (`src/core.c`) accepts exactly one argument (any other
arity returns `NULL`) and calls `activate_recursive`, which walks the expression tree
replacing every compound head `Inactive[h]` (an `Inactive` call of arity 1) by `h` and
rebuilding the node. The reactivated tree is handed back to the evaluator's fixed-point
loop, which then re-evaluates the now-active heads — so `Activate[Inactive[Integrate][g,
x]]` becomes `Integrate[g, x]` and actually integrates.

**Data structures.** Pure `Expr`-tree recursion. A non-function node is deep-copied
(`expr_copy`); a function node reactivates its head and every argument into a fresh
`Expr**` buffer that is passed to `expr_new_function` (which consumes the head and args).
The `Inactive[h] -> h` rewrite short-circuits before the generic head/argument walk, so
nested inert heads are peeled in one pass.

**Complexity / limits.** `O(size of expression)` for the single rewrite pass, plus
whatever re-evaluation the reactivated heads subsequently trigger. `Activate` is
`ATTR_PROTECTED` and is the exact inverse of `Inactive`.

- `Protected`.
- Reverses `Inactive`: `Activate[Inactive[Integrate][g, x]]` becomes `Integrate[g, x]` and evaluates.

**Attributes:** `Protected`.

## References

**See also:** [Inactive](../../expression-information/Inactive/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_inactive.c`](https://github.com/stblake/mathilda/blob/main/tests/test_inactive.c)

## Notes & additional examples

### Notes

`Activate[expr]` is the inverse of `Inactive`: it walks the expression replacing every
inert head `Inactive[h]` by `h` and then lets the evaluator re-run, so the reactivated
heads finally fire. It is the counterpart to holding a computation inert with `Inactive`.

The rewrite is recursive and reaches every nesting level, so a deeply inert tree is
reactivated in a single pass.
