# Context

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Context[] gives the current context ($Context).`**

**`Context[sym] gives the context in which sym resides.`**

**`Context["name"] gives the context of the symbol named "name" if it exists.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

```mathematica
In[1]:= Context[]
Out[1]= "Global`"

In[2]:= Context[Sin]
Out[2]= "System`"
```

## Implementation notes

`builtin_context` (`src/context.c`) has three forms. `Context[]` returns `context_current()` (the `g_current` string, default `"Global`"`). `Context[sym]` / `Context["name"]` (the symbol form is held — `ATTR_HOLDFIRST`) report the symbol's context: if the name already carries a backtick prefix it is peeled off and returned via `context_prefix_len`; otherwise the symbol is looked up with `symtab_lookup` and reported as `"System`"` when it is a builtin (`def->builtin_func != NULL`), else `"Global`"`. An unknown string-form name emits `Context::notfound` and leaves the call unevaluated; an unassigned bare symbol defaults to `"Global`"`.

**Attributes:** `HoldFirst`, `Protected`.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)
- Tests: [`tests/test_context.c`](https://github.com/stblake/mathilda/blob/main/tests/test_context.c)

## Notes & additional examples

### Notes

`Context[]` returns the current context (the same value as `$Context`), while `Context[sym]` returns the context in which a symbol lives — built-ins such as `Sin` reside in `"System`"`.
