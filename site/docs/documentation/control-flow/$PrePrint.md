# $PrePrint

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$PrePrint`**

is a global variable whose value, if set, is applied to every expression just before it is printed. Out\[n\] is assigned the unmodified result, but the printed form reflects the value returned by $PrePrint.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

```mathematica
In[1]:= $PrePrint = Framed
Out[1]= Framed[Framed]

In[2]:= 3 + 4
Out[2]= Framed[7]
```

## Implementation notes

A REPL session hook, not a builtin. Registered (docstring only) in `repl_hooks_init` (`src/repl_hooks.c`). `repl.c` calls `repl_apply_pre_print(out)` just before printing; if an OwnValue is set, `hook_call_eval` builds and evaluates `$PrePrint[expr]`. Crucially this is display-only: `Out[n]` is assigned the unmodified post-`$Post` result above, and only the rendered form reflects the `$PrePrint` value. Unset = identity.

**Attributes:** none registered.

## References

- Source: [`src/repl_hooks.c`](https://github.com/stblake/mathilda/blob/main/src/repl_hooks.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)

## Notes & additional examples

### Notes

`$PrePrint`, if set, is applied to every expression just before it is printed.
The displayed form reflects the `$PrePrint` value (here wrapped in `Framed`),
but `Out[n]` is assigned the unmodified result, so later references see the
plain value. Unset by default.
