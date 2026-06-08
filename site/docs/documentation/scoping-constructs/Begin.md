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
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= $Context
Out[1]= "Global`"

In[2]:= Begin["foo`"]
Out[2]= "foo`"

In[3]:= $Context
Out[3]= "foo`"

In[4]:= End[]
Out[4]= "foo`"

In[5]:= $Context
Out[5]= "Global`"
```

### Notes

`Begin["ctx`"]` switches `$Context` to `ctx`` and saves the previous context so the matching `End[]` can restore it. An argument that starts with a backtick is taken relative to the current context (e.g. `Begin["`Private`"]` inside `MyPkg`` yields `MyPkg`Private``).
