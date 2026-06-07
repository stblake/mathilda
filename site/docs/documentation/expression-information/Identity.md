# Identity

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Identity[expr] gives expr unchanged (the identity function).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Identity[x]
Out[1]= x

In[2]:= Identity[1 + 2]
Out[2]= 3

In[3]:= Map[Identity, {a, b, c}]
Out[3]= {a, b, c}
```

## Implementation notes

`builtin_identity` (`src/core.c`) is the one-argument identity: it returns a copy of its single argument unchanged, or `NULL` for any other arity.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
