# ValueQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ValueQ[expr]`**

gives True if a value has been defined for expr, False otherwise.

<details>
<summary>Notes</summary>

HoldAll: inspects the symbol itself, not its evaluated value. A bare symbol tests OwnValues; f\[...\] tests whether head f has any DownValues.

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= ValueQ[x]
Out[1]= False

In[2]:= x = 5; ValueQ[x]
Out[2]= True
```

Head has DownValues

```mathematica
In[3]:= f[x_] := x^2; {ValueQ[f[2]], ValueQ[f[a, b]]}
Out[3]= {True, True}
```

Bare symbol, only DownValues

```mathematica
In[4]:= ValueQ[f]
Out[4]= False
```

HoldAll preserved via Unevaluated

```mathematica
In[5]:= ValueQ /@ Unevaluated[{x, y}]
Out[5]= {True, False}
```

### Applications (4)

```mathematica
In[6]:= ValueQ[x]
Out[6]= False
```

Now x carries an OwnValue

```mathematica
In[7]:= x = 5; ValueQ[x]
Out[7]= True
```

True as soon as the head f has a DownValue

```mathematica
In[8]:= f[a_] := a^2; ValueQ[f[2]]
Out[8]= True
```

The bare symbol f has only DownValues, so no OwnValue

```mathematica
In[9]:= ValueQ[f]
Out[9]= False
```

## Implementation notes

**Algorithm.** `builtin_valueq` (`src/core.c`) inspects the symbol table rather
than evaluating its argument — `ValueQ` carries `HoldAll`, so the argument reaches
the builtin unevaluated. A bare `EXPR_SYMBOL` is valued iff its `SymbolDef`
carries a non-empty `own_values` chain (an immediate or delayed assignment). A
compound `f[...]` is valued iff the `SymbolDef` of its head `f` carries any
`down_values` — regardless of whether the particular arguments would match a rule.
Everything else (numbers, strings, undefined symbols) returns `False`.

**Data structures.** Two `symtab_lookup` probes at most, reading the `own_values`
/ `down_values` linked lists hanging off the `SymbolDef`; no rule matching and no
tree walk.

**Complexity / limits.** `O(1)` (a hash lookup plus a null test on a list head).
Because the test is existence-of-a-rule, not applicability, `ValueQ[f[a, b]]` is
`True` as soon as `f` has any `DownValue`. Wrong arity routes a `ValueQ::argx`
diagnostic through `builtin_arg_error`. Attributes `HoldAll`, `Protected`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_nminimize.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nminimize.c)

## Notes & additional examples

### Notes

`ValueQ[expr]` reports whether a value has been defined for `expr` **without
evaluating it** — it is `HoldAll`, so it inspects the symbol itself rather than the
value the symbol would evaluate to. A bare symbol is valued iff it carries an
`OwnValue` (an immediate `x = 5` or a delayed `y := RandomReal[]`); a compound
`f[...]` is valued iff its head `f` carries any `DownValue`, regardless of whether
the particular arguments match a rule — which is why `ValueQ[f[2]]` is `True` while
`ValueQ[f]` is `False`. Numbers, strings, and undefined symbols give `False`.
