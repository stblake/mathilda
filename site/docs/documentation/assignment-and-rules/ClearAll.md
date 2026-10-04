# ClearAll

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ClearAll[s1, s2, ...]`**

clears all values, definitions, attributes and messages for the named symbols. ClearAll\[{s1, s2, ...}\] accepts a list of specs.

<details>
<summary>Notes</summary>

ClearAll has attribute HoldAll; symbols with attribute Locked or Protected are not affected.

</details>

## Examples (7)

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

### Applications (4)

```mathematica
In[4]:= SetAttributes[g, Orderless]; g = 5
Out[4]= 5

In[5]:= ClearAll[g]

In[6]:= {g, Attributes[g]}
Out[6]= {g, {}}
```

Several symbols in one call

```mathematica
In[7]:= ClearAll[a, b]
```

## Implementation notes

**Algorithm.** `builtin_clear_all` (`src/core.c`) is one of the four
symbol-management heads driven by `core_apply_symbol_action`, which walks the
argument list and applies a per-symbol action to each spec — a bare symbol, a
string naming a symbol, or a flat `List` of such specs (so `ClearAll[{a, b}]`
works), with the name read out by `core_symbol_name_of`. The action here is
`core_clear_all_one`.

For each name, `core_clear_all_one` first skips any `Protected` or `Locked`
symbol (which is what shields every builtin), then does the full erase that
distinguishes `ClearAll` from `Clear`: `symtab_clear_symbol` drops the
OwnValues/DownValues, the attribute word is zeroed (bumping the rule epoch so
the evaluation cache is invalidated), and the docstring (usage message) is
freed. The C builtin function pointer, if any, is left intact — the `Protected`
guard already keeps `ClearAll` away from builtins. The head returns `Null`.

**Attributes & limits.** `ClearAll` carries `HoldAll | Protected`, so the
symbols arrive unevaluated rather than being replaced by their current values.
Non-symbol/non-string specs are silently ignored.

- `ClearAll` has attributes `{HoldAll, Protected}`; `Remove` has
  `{HoldAll, Locked, Protected}`. Both hold their arguments, so they operate on
  the symbol, not its current value.
- Neither affects symbols with the attribute `Locked` or `Protected`. This is
  what prevents `Remove`/`ClearAll` from ever deleting or wiping a built-in.
- `ClearAll`, unlike `Clear`, also removes attributes and the usage message.
- Both return `Null`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Remove](../../assignment-and-rules/Remove/), [Clear](../../assignment-and-rules/Clear/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)
- Tests: [`tests/test_condition_downvalue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_condition_downvalue.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_evaluate.c`](https://github.com/stblake/mathilda/blob/main/tests/test_evaluate.c)

## Notes & additional examples

### Notes

`ClearAll[s]` is the thorough erase: it removes `s`'s OwnValues and DownValues
*and* its attributes and usage message, so the symbol is returned to its
pristine undefined state — compare `Clear[s]`, which drops only the
values. Here `g` loses both its assigned value and the `Orderless` attribute, so
`Attributes[g]` is `{}`.

Arguments may be symbols, strings naming symbols, or a flat list of them
(`ClearAll[{a, b}]`). A `Protected` or `Locked` symbol is skipped, which is what
keeps `ClearAll` from ever gutting a built-in. The result is `Null`.
