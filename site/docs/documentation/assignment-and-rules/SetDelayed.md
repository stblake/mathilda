# SetDelayed

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
lhs := rhs or SetDelayed[lhs, rhs]
    assigns rhs to lhs as a delayed rule: rhs is held and evaluated
    each time the rule fires (with bindings from lhs substituted in),
    not at assignment time.  SetDelayed has attribute HoldAll.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `SetDelayed` (`:=`) shares the special-primitive path with `Set` in `evaluate_step` (`src/eval.c`, step 6, `head_name == SYM_SetDelayed`). The head carries `ATTR_HOLDALL` (`src/attr.c`), so neither side is pre-evaluated: the RHS is stored *unevaluated* and re-evaluated each time the rule fires. As with `Set`, the LHS is partially evaluated to canonicalise the assignment target (arguments evaluated, head and pattern-bearing/destructuring sub-elements held — see `lhs_arg_contains_pattern`), then `apply_assignment(lhs, rhs, is_delayed=true)` installs the value. `SetDelayed` returns the symbol `Null` rather than the RHS.

**Data structures.** Symbol LHS → OwnValue (`symtab_add_own_value`); function-head LHS → DownValue (`symtab_add_down_value`) keyed on the head, with the full LHS as the match pattern and the held RHS as the replacement; `List` LHS recurses for destructuring. One delayed-specific transform: if the RHS is `Condition[body, test]`, `apply_assignment` rewrites the rule as `Condition[lhs, test] := body`, so `f[x_] := body /; test` is stored equivalently to `f[x_] /; test := body`; the condition's head symbol is then unwrapped to find the true DownValue key. Protected targets are refused with a `SetDelayed::wrsym` message. For `$RecursionLimit` a copy of the RHS is evaluated to sync the C-side limit.

**Attributes:** `HoldAll`, `Protected`, `SequenceHold`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= f[n_] := n^2
Out[1]= Null

In[2]:= f[4]
Out[2]= 16

In[3]:= f[a + b]
Out[3]= (a + b)^2
```

### Notes

`:=` (`SetDelayed`) holds the right-hand side and re-evaluates it each time the rule fires, with the pattern bindings substituted in. The definition itself returns `Null`.
