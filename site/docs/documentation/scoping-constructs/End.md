# End

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`End[] restores the context that was active before the matching Begin[]`**

<details>
<summary>Notes</summary>

and returns the closed context as a string.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= $Context
Out[1]= "Global`"

In[2]:= Begin["foo`"]
Out[2]= "foo`"

In[3]:= End[]
Out[3]= "foo`"

In[4]:= $Context
Out[4]= "Global`"
```

## Implementation notes

`builtin_end` (0-arg) calls `context_end`, which pops the top frame of the context stack (`g_stack`, a linked list of `CtxFrame` pushed by `Begin`/`BeginPackage`). `frame_pop` restores `$Context` (`g_current`) and the `$ContextPath` snapshot (`saved_path`) that the matching `Begin[]` saved, frees the popped frame, and republishes the live state. The builtin returns the *closed* context string (the one that was current just before the pop). If the stack is empty it emits `End::noctx` and returns NULL (unevaluated). Context names are owned as plain `char*` via the file's `ctx_strdup` (a C99-safe `strdup` replacement).

**Attributes:** `Protected`.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)
- Tests: [`tests/test_context.c`](https://github.com/stblake/mathilda/blob/main/tests/test_context.c)
- Tests: [`tests/test_symbol.c`](https://github.com/stblake/mathilda/blob/main/tests/test_symbol.c)

## Notes & additional examples

### Notes

`End[]` closes the context opened by the matching `Begin[]`, restoring the previously active `$Context` and returning the context it just closed as a string.
