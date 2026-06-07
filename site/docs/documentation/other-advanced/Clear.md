# Clear

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Clear[s1, s2, ...]
    clears all OwnValues and DownValues attached to the named symbols,
    leaving attributes and the symbol itself intact.
Clear has attribute HoldAll; Protected symbols are skipped with a
diagnostic.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_clear` (`src/core.c`) iterates its arguments and, for each that is a symbol, calls `symtab_clear_symbol(name)` to remove that symbol's OwnValues and DownValues (its rules/assignments) while leaving the symbol itself, its attributes, and any builtin binding intact. Non-symbol arguments are ignored. Returns `Null`. It carries `ATTR_HOLDALL` so the symbols are not evaluated to their current values before being cleared.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= x = 5
Out[1]= 5

In[2]:= x + 1
Out[2]= 6

In[3]:= Clear[x]
Out[3]= Null

In[4]:= x + 1
Out[4]= 1 + x
```

### Notes

`Clear[s]` removes all OwnValues and DownValues attached to `s`, so the symbol becomes undefined again (attributes and the symbol itself are left intact).
