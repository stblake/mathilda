# Remove

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Remove[s1, s2, ...]`**

removes the named symbols completely, deleting their definitions from the symbol table. Remove\[{s1, s2, ...}\] accepts a list of specs.

<details>
<summary>Notes</summary>

Remove has attribute HoldAll; symbols with attribute Locked or Protected are not affected.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

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
  `{HoldAll, Locked, Protected}`. Both hold their arguments, so they operate on
  the symbol, not its current value.
- Neither affects symbols with the attribute `Locked` or `Protected`. This is
  what prevents `Remove`/`ClearAll` from ever deleting or wiping a built-in.
- `ClearAll`, unlike `Clear`, also removes attributes and the usage message.
- Both return `Null`.

**Attributes:** `HoldAll`, `Locked`, `Protected`.

## References

**See also:** [ClearAll](../../assignment-and-rules/ClearAll/), [Clear](../../assignment-and-rules/Clear/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)
