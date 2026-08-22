# BeginPackage

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BeginPackage["ctx`"] sets the current context to "ctx`" and restricts`**

**`BeginPackage["ctx`", {"need1`", ...}] additionally prepends the`**

<details>
<summary>Notes</summary>

$ContextPath to {"ctx\`", "System\`"}, matching the standard package prologue. listed contexts to $ContextPath.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= $Context
Out[1]= "Global`"

In[2]:= BeginPackage["MyPkg`"]
Out[2]= "MyPkg`"

In[3]:= $Context
Out[3]= "MyPkg`"

In[4]:= $ContextPath
Out[4]= {"MyPkg`", "System`"}
```

## Implementation notes

`builtin_begin_package` (`src/context.c`) accepts `"ctx`"` and an optional `Needs` list of context strings, dispatching to `context_begin_package`. That pushes a `FRAME_PACKAGE` frame (snapshotting context and path), sets `g_current` to the absolute package context (relative `` `... `` contexts are rejected), clears the search path and rebuilds it as `{ctx`, System`}` plus any valid, non-duplicate `Needs` entries. `republish_state` then republishes `$Context`/`$ContextPath`, and the new context is returned. `EndPackage[]` (`context_end_package`) pops the frame and prepends the closed context to `$ContextPath`. Invalid specs emit `BeginPackage::cxt`.

**Attributes:** `Protected`.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)
- Tests: [`tests/test_context.c`](https://github.com/stblake/mathilda/blob/main/tests/test_context.c)

## Notes & additional examples

### Notes

`BeginPackage["ctx`"]` sets `$Context` to `ctx`` and restricts `$ContextPath` to `{"ctx`", "System`"}`, matching the standard package prologue so only system symbols and the package's own symbols resolve under short names. Use the matching `EndPackage[]` to close it.
