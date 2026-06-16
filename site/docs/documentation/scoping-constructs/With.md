# With

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
With[{x = x0, ...}, expr] specifies that x should be replaced by x0 throughout expr.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= x = 10; With[{x = 5}, x^2]
Out[1]= 25
```

## Implementation notes

**Algorithm.** `builtin_with` (in `src/modular.c`) implements lexical constants by **literal substitution** — no renaming and no symbol-table mutation (contrast `Module`/`Block`). `With` carries `HoldAll | Protected` (set in `src/attr.c`). Each binding must be `x = val` (value evaluated immediately in the outer scope) or `x := val` (RHS substituted verbatim, unevaluated); the handler collects these into a `ScopingEnv` linked list mapping the constant name to its replacement expression.

The body is rewritten by `substitute_scoping`, the same recursive, shadow-aware tree walk used by `Module`: it replaces free occurrences of each constant, drops a name from the environment when descending into a nested scoping construct that rebinds it, and substitutes into nested binding RHSs (so `With[{q=12}, With[{k=q}, k]]` resolves `k` to 12) without touching binding LHS names. The substituted body is then `evaluate`d; `Return[v]`/`Return[v, With]` is trapped via `eval_classify_return`. Because substitution is structural, the constants vanish before evaluation and leave no trace in the symbol table.

- `HoldAll`, `Protected`.
- Replaces occurrences of symbols in the body before evaluation.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/modular.c`](https://github.com/stblake/mathilda/blob/main/src/modular.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= With[{a = 2, b = 3}, a^2 + b^2]
Out[1]= 13
```

The bound values are substituted *before* the body evaluates, so injecting an exact algebraic constant lets the symbolic engine verify a defining identity — here the golden ratio satisfying `phi^2 - phi - 1 == 0`:

```mathematica
In[1]:= With[{phi = (1 + Sqrt[5])/2}, Simplify[phi^2 - phi - 1]]
Out[1]= 0
```

`With` localizes the value cleanly for numeric work too, fixing a parameter and driving a high-precision computation:

```mathematica
In[1]:= With[{n = 20}, N[Sum[1/k^2, {k, 1, n}], 30]]
Out[1]= 1.596163243913023316640878872058
```

Bindings flow through complex arithmetic transparently:

```mathematica
In[1]:= With[{x = 1 + I}, x^2]
Out[1]= 2*I
```

### Notes

`With[{x = x0, ...}, expr]` replaces each `x` by its value `x0` throughout `expr` and then evaluates the result. Unlike `Module`, the substitution is literal and immediate, making `With` ideal for inlining constants and parameters.
