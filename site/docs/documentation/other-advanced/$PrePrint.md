# $PrePrint

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$PrePrint
    is a global variable whose value, if set, is applied to every
    expression just before it is printed. Out[n] is assigned the
    unmodified result, but the printed form reflects the value
    returned by $PrePrint.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

A REPL session hook, not a builtin. Registered (docstring only) in `repl_hooks_init` (`src/repl_hooks.c`). `repl.c` calls `repl_apply_pre_print(out)` just before printing; if an OwnValue is set, `hook_call_eval` builds and evaluates `$PrePrint[expr]`. Crucially this is display-only: `Out[n]` is assigned the unmodified post-`$Post` result above, and only the rendered form reflects the `$PrePrint` value. Unset = identity.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/repl_hooks.c`](https://github.com/stblake/mathilda/blob/main/src/repl_hooks.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
