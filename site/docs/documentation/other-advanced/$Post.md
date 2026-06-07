# $Post

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$Post
    is a global variable whose value, if set, is applied to every
    output expression after evaluation.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

A REPL session hook, not a builtin. Registered (docstring only, no default OwnValue) in `repl_hooks_init` (`src/repl_hooks.c`). In each REPL cycle `repl.c` calls `repl_apply_post(result)` on the evaluator's output; when an OwnValue is set the helper builds `$Post[expr]` and runs it through `evaluate()` via `hook_call_eval`. This happens after evaluation and before `Out[n]` is stored, so `$Post` can transform the visible result. Unset = identity.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/repl_hooks.c`](https://github.com/stblake/mathilda/blob/main/src/repl_hooks.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= $Post = Framed
Out[1]= Framed[Framed]

In[2]:= 2 + 2
Out[2]= Framed[4]
```

### Notes

`$Post`, if set, is applied to every output expression after evaluation. Here
`$Post = Framed` wraps each result in `Framed[...]` (note the assignment's own
echo is also wrapped). Unset by default, in which case results are shown
unmodified.
