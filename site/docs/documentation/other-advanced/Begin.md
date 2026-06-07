# Begin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Begin["ctx`"] sets the current context ($Context) to "ctx`", saving
the previous value for End[] to restore. If the argument starts with a
backtick it is interpreted relative to the current context:
Begin["`Private`"] inside "MyPkg`" yields "MyPkg`Private`".
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_begin` (`src/context.c`) takes a single string `"ctx`"` and calls `context_begin`, which pushes a `FRAME_BEGIN` onto the internal `g_stack` (snapshotting the current context and search path), then sets `g_current` to the new context. A leading-backtick argument like `"`Private`"` is interpreted relative to the current context (`MyPkg``` + `Private``). `republish_state` refreshes the `$Context`/`$ContextPath` OwnValues, and the new current context string is returned. Invalid (non-backtick-terminated) specs emit `Begin::cxt`. The matching `End[]` (`builtin_end`) pops the frame to restore the saved context.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
