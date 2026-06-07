# EndPackage

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
EndPackage[] restores the state saved by BeginPackage and prepends the
just-closed package context to $ContextPath so its exported symbols are
visible under short names. Returns Null.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_end_package` (0-arg) calls `context_end_package`, which pops the top `CtxFrame` from the context stack via `frame_pop` (restoring `$Context` and the saved `$ContextPath`), then *prepends* the just-closed package context to `$ContextPath` with `path_prepend` (unless already present) so symbols defined inside the package remain visible after the package closes. State is republished. Unlike `End[]`, the builtin returns the symbol `Null`. On an empty stack it emits `EndPackage::noctx` and returns NULL. The frame stack is shared by `Begin`/`BeginPackage`/`End`/`EndPackage`; the `FrameKind` tag distinguishes package frames but pop logic is common.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= BeginPackage["MyPkg`"]
Out[1]= "MyPkg`"

In[2]:= EndPackage[]
Out[2]= Null

In[3]:= $Context
Out[3]= "Global`"

In[4]:= $ContextPath
Out[4]= {"MyPkg`", "Global`", "System`"}
```

### Notes

`EndPackage[]` restores the context state saved by `BeginPackage` and prepends the just-closed package context to `$ContextPath`, so the package's exported symbols remain visible under their short names. It returns `Null`.
