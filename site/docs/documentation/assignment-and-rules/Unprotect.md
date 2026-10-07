# Unprotect

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Unprotect[s1, s2, ...]`**

removes the attribute Protected from the named symbols and returns the list of their names. Unprotect\[{s1, ...}\] accepts a list of specs.

<details>
<summary>Notes</summary>

Unprotect has attribute HoldAll.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= f[x_] := x^2; Protect[f]
Out[1]= {"f"}

In[2]:= Unprotect[f]
Out[2]= {"f"}
```

### Applications (3)

```mathematica
In[3]:= Protect[baz]
Out[3]= {"baz"}
```

Returns the names whose Protected attribute was actually cleared

```mathematica
In[4]:= Unprotect[baz]
Out[4]= {"baz"}
```

```mathematica
In[5]:= MemberQ[Attributes[baz], Protected]
Out[5]= False
```

## Implementation notes

**Algorithm.** `builtin_unprotect` (`src/core.c`) is the mirror of `Protect`: it
calls the shared driver `core_protect_unprotect(res, protecting = false)`, which
walks the argument list (symbols, strings, or a flat `List` of them) and applies
`core_unprotect_one` to each name. `core_unprotect_one` does nothing if the
symbol is not currently `Protected`; otherwise it
clears the `ATTR_PROTECTED` bit and bumps the rule epoch, returning `true` only
when the bit was actually cleared.

The driver returns a `List` of the names (as strings) whose protection state
changed, so `Unprotect` of a symbol that was never protected yields `{}`. Once
the bit is cleared, `Set`/`SetDelayed` and the clearing heads will again accept
the symbol, which is the usual prelude to redefining or extending a built-in's
behaviour.

**Attributes & limits.** `Unprotect` carries `HoldAll | Protected`, so its
arguments arrive as unevaluated names.

- Both have attributes `{HoldAll, Protected}` and hold their arguments.
- The typical sequence for adding rules to an existing symbol is
  `Unprotect[f]; definition; Protect[f]`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Protect](../../assignment-and-rules/Protect/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)

## Notes & additional examples

### Notes

`Unprotect[s1, s2, ...]` is the inverse of `Protect`: it clears the `Protected`
attribute from each symbol and returns the list of names (as strings) whose state
changed, so unprotecting a symbol that was never protected gives `{}`.

Once a symbol is unprotected, `Set`/`SetDelayed` and the clearing heads accept it
again — the standard way to override or extend a built-in's definition.
Arguments may be symbols, strings, or a flat list.
