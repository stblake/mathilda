# Unprotect

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Unprotect[s1, s2, ...]
    removes the attribute Protected from the named symbols and returns
    the list of their names. Unprotect[{s1, ...}] accepts a list of specs.
Unprotect has attribute HoldAll; Locked symbols are not affected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= f[x_] := x^2; Protect[f]
Out[1]= {"f"}

In[2]:= Unprotect[f]
Out[2]= {"f"}
```

## Implementation notes

- Both have attributes `{HoldAll, Protected}` and hold their arguments.
- Neither affects symbols with the attribute `Locked`.
- The typical sequence for adding rules to an existing symbol is

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
