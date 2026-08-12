# Clear

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Clear[s1, s2, ...]`**

clears all OwnValues and DownValues attached to the named symbols, leaving attributes and the symbol itself intact.

<details>
<summary>Notes</summary>

Clear has attribute HoldAll; Protected symbols are skipped with a diagnostic.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

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

## Implementation notes

`builtin_clear` (`src/core.c`) iterates its arguments and, for each that is a symbol, calls `symtab_clear_symbol(name)` to remove that symbol's OwnValues and DownValues (its rules/assignments) while leaving the symbol itself, its attributes, and any builtin binding intact. Non-symbol arguments are ignored. Returns `Null`. It carries `ATTR_HOLDALL` so the symbols are not evaluated to their current values before being cleared.

**Attributes:** `HoldAll`, `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_catch_throw.c`](https://github.com/stblake/mathilda/blob/main/tests/test_catch_throw.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
- Tests: [`tests/test_cond.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cond.c)

## Notes & additional examples

### Notes

`Clear[s]` removes all OwnValues and DownValues attached to `s`, so the symbol becomes undefined again (attributes and the symbol itself are left intact).
