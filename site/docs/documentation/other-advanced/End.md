# End

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
End[] restores the context that was active before the matching Begin[]
and returns the closed context as a string.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_end` (0-arg) calls `context_end`, which pops the top frame of the context stack (`g_stack`, a linked list of `CtxFrame` pushed by `Begin`/`BeginPackage`). `frame_pop` restores `$Context` (`g_current`) and the `$ContextPath` snapshot (`saved_path`) that the matching `Begin[]` saved, frees the popped frame, and republishes the live state. The builtin returns the *closed* context string (the one that was current just before the pop). If the stack is empty it emits `End::noctx` and returns NULL (unevaluated). Context names are owned as plain `char*` via the file's `ctx_strdup` (a C99-safe `strdup` replacement).

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
