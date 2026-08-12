# EndPackage

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EndPackage[] restores the state saved by BeginPackage and prepends the`**

<details>
<summary>Notes</summary>

just-closed package context to $ContextPath so its exported symbols are visible under short names. Returns Null.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

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

## Implementation notes

`builtin_end_package` (0-arg) calls `context_end_package`, which pops the top `CtxFrame` from the context stack via `frame_pop` (restoring `$Context` and the saved `$ContextPath`), then *prepends* the just-closed package context to `$ContextPath` with `path_prepend` (unless already present) so symbols defined inside the package remain visible after the package closes. State is republished. Unlike `End[]`, the builtin returns the symbol `Null`. On an empty stack it emits `EndPackage::noctx` and returns NULL. The frame stack is shared by `Begin`/`BeginPackage`/`End`/`EndPackage`; the `FrameKind` tag distinguishes package frames but pop logic is common.

**Attributes:** `Protected`.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)
- Tests: [`tests/test_context.c`](https://github.com/stblake/mathilda/blob/main/tests/test_context.c)

## Notes & additional examples

### Notes

`EndPackage[]` restores the context state saved by `BeginPackage` and prepends the just-closed package context to `$ContextPath`, so the package's exported symbols remain visible under their short names. It returns `Null`.
