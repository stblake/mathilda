# Unprotect

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Unprotect[s1, s2, ...]`**

removes the attribute Protected from the named symbols and returns the list of their names. Unprotect\[{s1, ...}\] accepts a list of specs.

<details>
<summary>Notes</summary>

Unprotect has attribute HoldAll; Locked symbols are not affected.

</details>

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

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
  `Unprotect[f]; definition; Protect[f]`.

**Attributes:** `HoldAll`, `Protected`.

## See also

[Protect](../../assignment-and-rules/Protect/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)
