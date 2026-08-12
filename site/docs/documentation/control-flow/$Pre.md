# $Pre

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$Pre`**

is a global variable whose value, if set, is applied to every input expression after parsing and before standard evaluation.

<details>
<summary>Notes</summary>

Unless $Pre is assigned to a head with HoldAll, the wrapped expression is evaluated before $Pre sees it -- which makes the effect indistinguishable from $Post.

</details>

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

```mathematica
In[1]:= $Pre = Hold
Out[1]= Hold

In[2]:= 1 + 1
Out[2]= Hold[1 + 1]
```

## Implementation notes

A REPL session hook, not a builtin. `repl_hooks_init` (`src/repl_hooks.c`) merely touches the symbol so `?$Pre` works; no default OwnValue is installed, so out of the box the hook is a no-op. Each REPL iteration `repl.c` calls `repl_apply_pre(parsed)`; if `hook_is_set("$Pre")` (i.e. `symtab_get_own_values("$Pre")` is non-empty) it builds `$Pre[expr]` and runs it through the standard `evaluate()` via `hook_call_eval`, applied after parsing but before the main evaluation. Because the wrapped expression is evaluated before `$Pre` sees it unless `$Pre` is assigned a `HoldAll` head, its effect is usually indistinguishable from `$Post`.

**Attributes:** none registered.

## References

- Source: [`src/repl_hooks.c`](https://github.com/stblake/mathilda/blob/main/src/repl_hooks.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)

## Notes & additional examples

### Notes

`$Pre`, if set, is applied to every input expression after parsing and before
standard evaluation. Using a head with `HoldAll` (such as `Hold`) lets `$Pre`
intercept the unevaluated input; otherwise the argument is evaluated first and
the effect is indistinguishable from `$Post`. Unset by default.
