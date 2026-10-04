# SymbolName

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`SymbolName[symbol] gives the name of symbol as a string, with any context prefix removed.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

The symbol's own name as a string

```mathematica
In[1]:= SymbolName[xyz]
Out[1]= "xyz"
```

The context prefix is stripped

```mathematica
In[2]:= SymbolName[Global`foo]
Out[2]= "foo"
```

The result is a String

```mathematica
In[3]:= Head[SymbolName[abc]]
Out[3]= String
```

## Implementation notes

**Algorithm.** `builtin_symbolname` takes one argument and, when it is a symbol,
returns the symbol's short name as a string with any context prefix stripped. The
implementation finds the last backtick in the interned name with `strrchr(n, '`')`
and returns everything after it (or the whole name if there is no backtick), so
`SymbolName[Global`x]` and `SymbolName[P`Private`x]` both give `"x"`. Anything that
is not a symbol — a number, a compound expression — declines (returns `NULL`) and
the call is left unevaluated, which is the same observable a caller testing the
result's head would get from Wolfram's `SymbolName::sym` path.

**Data structures.** None beyond the argument. The returned value is a freshly
allocated `EXPR_STRING` copied from the tail of the interned name; the name buffer
itself is not modified.

**Complexity / limits.** `O(k)` in the length of the name (one backtick scan). The
argument is evaluated first, so `SymbolName[s]` where `s` has a value reports on
that value, not on `s`; apply it to the bare symbol to name the symbol itself.
Wrong arity declines.

**Attributes:** `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`SymbolName[sym]` gives the short name of `sym` as a string, with any context
prefix (`Global\``, `P\`Private\``, ...) removed — only the part after the last
backtick is kept. The argument is evaluated first, so to name a symbol that has a
value, pass it held or use the bare symbol. A non-symbol argument is left
unevaluated.
