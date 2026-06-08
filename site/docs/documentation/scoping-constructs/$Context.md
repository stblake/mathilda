# $Context

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$Context is a string giving the current context. New symbols are created
in this context unless otherwise qualified. Modified by Begin[],
BeginPackage[], End[], and EndPackage[].
```

## Examples

_No verified examples yet for this function._

## Implementation notes

Not a builtin: `$Context` is an ordinary symbol whose OwnValue mirrors the internal context state held in `src/context.c` (`g_current`, a heap string always ending in a backtick). Whenever `Begin`/`BeginPackage`/`End`/`EndPackage` change state, `republish_state` calls `publish_own("$Context", context_as_string())`, which clears the symbol's OwnValues and reinstalls a single string-valued rule. The symbol is marked `ATTR_PROTECTED` in `context_init` so the user cannot reassign it directly — mutation only flows through the context-stack builtins. `context_resolve_name` reads `g_current` to qualify bare identifiers at parse/lookup time.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= $Context
Out[1]= "Global`"

In[2]:= Begin["foo`"]
Out[2]= "foo`"

In[3]:= $Context
Out[3]= "foo`"
```

### Notes

`$Context` is the string giving the current context; new unqualified symbols are created here. It is read-only for direct assignment — change it through `Begin`, `BeginPackage`, `End`, and `EndPackage`.
