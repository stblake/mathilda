# Remove

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Remove[s1, s2, ...]
    removes the named symbols completely, deleting their definitions
    from the symbol table. Remove[{s1, s2, ...}] accepts a list of specs.
Remove has attribute HoldAll; symbols with attribute Locked or
Protected are not affected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= f[x_] := x^2; SetAttributes[f, Listable]; Attributes[f]
Out[1]= {Listable}

In[2]:= ClearAll[f]; {Attributes[f], DownValues[f]}
Out[2]= {{}, {}}

In[3]:= x = 2; Remove[x]; x
Out[3]= x
```

## Implementation notes

- `ClearAll` has attributes `{HoldAll, Protected}`; `Remove` has

**Attributes:** `HoldAll`, `Locked`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
