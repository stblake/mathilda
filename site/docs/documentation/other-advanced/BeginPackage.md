# BeginPackage

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
BeginPackage["ctx`"] sets the current context to "ctx`" and restricts
$ContextPath to {"ctx`", "System`"}, matching Mathematica's package
prologue.
BeginPackage["ctx`", {"need1`", ...}] additionally prepends the
listed contexts to $ContextPath.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_begin_package` (`src/context.c`) accepts `"ctx`"` and an optional `Needs` list of context strings, dispatching to `context_begin_package`. That pushes a `FRAME_PACKAGE` frame (snapshotting context and path), sets `g_current` to the absolute package context (relative `` `... `` contexts are rejected), clears the search path and rebuilds it as `{ctx`, System`}` plus any valid, non-duplicate `Needs` entries. `republish_state` then republishes `$Context`/`$ContextPath`, and the new context is returned. `EndPackage[]` (`context_end_package`) pops the frame and prepends the closed context to `$ContextPath`. Invalid specs emit `BeginPackage::cxt`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

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

### Notes

`BeginPackage["ctx`"]` sets `$Context` to `ctx`` and restricts `$ContextPath` to `{"ctx`", "System`"}`, matching Mathematica's package prologue so only system symbols and the package's own symbols resolve under short names. Use the matching `EndPackage[]` to close it.
