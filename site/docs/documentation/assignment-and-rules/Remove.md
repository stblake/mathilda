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

## Examples (6)

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

### Applications (3)

```mathematica
In[4]:= tmpvar = 7
Out[4]= 7

In[5]:= Remove[tmpvar]

In[6]:= tmpvar
Out[6]= tmpvar
```

## Implementation notes

**Algorithm.** `builtin_remove` (`src/core.c`) shares the `core_apply_symbol_action`
walker with `ClearAll`: it visits each argument — a symbol, a string, or a flat
`List` of them — and applies `core_remove_one` to every resolved name, returning
`Null`. `core_remove_one` skips any `Protected` or `Locked` symbol (the guard
that keeps `Remove` from ever deleting a builtin) and otherwise calls
`symtab_remove_symbol(name)`, deleting the symbol's definition from the symbol
table entirely.

This is a stronger erase than `ClearAll`: where `ClearAll` empties a symbol but
leaves the entry in place, `Remove` deletes the entry, so the name no longer
appears in the symbol table until it is next referenced (at which point a fresh,
undefined symbol is created). `Remove` is itself registered `Locked` so it cannot
be removed.

**Attributes & limits.** `Remove` carries `HoldAll | Locked | Protected`; its
arguments arrive unevaluated, and non-symbol/non-string specs are ignored.

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

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)

## Notes & additional examples

### Notes

`Remove[s]` deletes the symbol `s` from the symbol table entirely — a stronger
erase than `ClearAll`, which empties a symbol but keeps its entry. After removal
the name no longer exists; the next reference to it (as in `In[3]`) creates a
fresh, undefined symbol. The result is `Null`.

Arguments may be symbols, strings, or a flat list of them. A `Protected` or
`Locked` symbol is skipped, so `Remove` can never delete a built-in; `Remove`
itself is `Locked`.
