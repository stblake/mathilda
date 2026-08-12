# Begin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Begin["ctx`"] sets the current context ($Context) to "ctx`", saving`**

**`Begin["`Private`"] inside "MyPkg`" yields "MyPkg`Private`".`**

<details>
<summary>Notes</summary>

the previous value for End\[\] to restore. If the argument starts with a backtick it is interpreted relative to the current context:

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

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

## Implementation notes

`builtin_begin` (`src/context.c`) takes a single string `"ctx`"` and calls `context_begin`, which pushes a `FRAME_BEGIN` onto the internal `g_stack` (snapshotting the current context and search path), then sets `g_current` to the new context. A leading-backtick argument like `"`Private`"` is interpreted relative to the current context (`MyPkg``` + `Private``). `republish_state` refreshes the `$Context`/`$ContextPath` OwnValues, and the new current context string is returned. Invalid (non-backtick-terminated) specs emit `Begin::cxt`. The matching `End[]` (`builtin_end`) pops the frame to restore the saved context.

**Attributes:** `Protected`.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)
- Tests: [`tests/test_context.c`](https://github.com/stblake/mathilda/blob/main/tests/test_context.c)
- Tests: [`tests/test_symbol.c`](https://github.com/stblake/mathilda/blob/main/tests/test_symbol.c)

## Notes & additional examples

### Notes

`Begin["ctx`"]` switches `$Context` to `ctx`` and saves the previous context so the matching `End[]` can restore it. An argument that starts with a backtick is taken relative to the current context (e.g. `Begin["`Private`"]` inside `MyPkg`` yields `MyPkg`Private``).
