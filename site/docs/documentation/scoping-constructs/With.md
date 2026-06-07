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
