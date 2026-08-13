# $ContextPath

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

$ContextPath is a list of contexts used (in order) to resolve bare identifiers to existing qualified symbols. Modified by BeginPackage\[\] and EndPackage\[\].

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= $ContextPath
Out[1]= {"Global`", "System`"}

In[2]:= BeginPackage["MyPkg`"]
Out[2]= "MyPkg`"

In[3]:= $ContextPath
Out[3]= {"MyPkg`", "System`"}
```

## Implementation notes

Not a builtin: `$ContextPath` is a symbol whose OwnValue tracks the internal search path `g_path` (a `char**` of context strings, most-specific first) in `src/context.c`. `republish_state` rebuilds it via `publish_own("$ContextPath", context_path_as_list())`, which materialises the `g_path` entries into a `List` of strings. It is `ATTR_PROTECTED`; the list is only mutated by `BeginPackage` (which resets the path to `{ctx`, System`}` plus any `Needs` contexts) and `EndPackage` (which prepends the just-closed package context). `context_resolve_name` walks this path in order to resolve a bare name to an existing qualified symbol.

**Attributes:** `Protected`.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)

## Notes & additional examples

### Notes

`$ContextPath` is the ordered list of contexts searched to resolve bare identifiers to existing qualified symbols. It is read-only for direct assignment and is modified by `BeginPackage` and `EndPackage`.
