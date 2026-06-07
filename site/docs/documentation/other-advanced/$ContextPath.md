# $ContextPath

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$ContextPath is a list of contexts used (in order) to resolve bare
identifiers to existing qualified symbols. Modified by BeginPackage[]
and EndPackage[].
```

## Examples

_No verified examples yet for this function._

## Implementation notes

Not a builtin: `$ContextPath` is a symbol whose OwnValue tracks the internal search path `g_path` (a `char**` of context strings, most-specific first) in `src/context.c`. `republish_state` rebuilds it via `publish_own("$ContextPath", context_path_as_list())`, which materialises the `g_path` entries into a `List` of strings. It is `ATTR_PROTECTED`; the list is only mutated by `BeginPackage` (which resets the path to `{ctx`, System`}` plus any `Needs` contexts) and `EndPackage` (which prepends the just-closed package context). `context_resolve_name` walks this path in order to resolve a bare name to an existing qualified symbol.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
