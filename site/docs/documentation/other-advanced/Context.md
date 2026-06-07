# Context

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Context[] gives the current context ($Context).
Context[sym] gives the context in which sym resides.
Context["name"] gives the context of the symbol named "name" if it exists.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_context` (`src/context.c`) has three forms. `Context[]` returns `context_current()` (the `g_current` string, default `"Global`"`). `Context[sym]` / `Context["name"]` (the symbol form is held — `ATTR_HOLDFIRST`) report the symbol's context: if the name already carries a backtick prefix it is peeled off and returned via `context_prefix_len`; otherwise the symbol is looked up with `symtab_lookup` and reported as `"System`"` when it is a builtin (`def->builtin_func != NULL`), else `"Global`"`. An unknown string-form name emits `Context::notfound` and leaves the call unevaluated; an unassigned bare symbol defaults to `"Global`"`.

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
