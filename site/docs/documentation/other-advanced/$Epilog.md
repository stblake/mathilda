# $Epilog

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$Epilog
    is a symbol whose value, if any, is evaluated once when the
    Mathilda session terminates (via Quit[] or EOF).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

A REPL session hook, not a builtin. Registered (docstring only) in `repl_hooks_init` (`src/repl_hooks.c`). Unlike the per-line hooks, `$Epilog` is a bare symbol evaluated once at session teardown: `repl.c` calls `repl_apply_epilog()` on `Quit[]` or EOF, which — if an OwnValue is set — evaluates the symbol `$Epilog` (not a call) for its side effects and discards the result. Unset = no-op.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/repl_hooks.c`](https://github.com/stblake/mathilda/blob/main/src/repl_hooks.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= $Epilog := Print["bye"]
Out[1]= Null

In[2]:= Quit[]
"bye"
```

### Notes

`$Epilog` is a symbol whose value, if any, is evaluated exactly once when the
session terminates via `Quit[]` or EOF. Assigning it with `:=` (delayed) lets
it run cleanup or a parting action at shutdown; the `"bye"` above is printed as
the session exits.
