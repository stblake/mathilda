# BooleanQ

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`BooleanQ[expr] gives True if expr is either True or False, and False otherwise. Unlike TrueQ it tests the symbol, so BooleanQ[False] is True.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= BooleanQ[True]
Out[1]= True
```

The inequality reduces to True first

```mathematica
In[2]:= BooleanQ[2 > 1]
Out[2]= True
```

A plain symbol is not a boolean

```mathematica
In[3]:= BooleanQ[x]
Out[3]= False
```

Nor is a number

```mathematica
In[4]:= BooleanQ[1]
Out[4]= False
```

## Implementation notes

**Algorithm.** `builtin_booleanq` is a one-argument predicate. It checks that its
single argument is an `EXPR_SYMBOL` whose interned name is either `SYM_True` or
`SYM_False`, and returns the symbol `True` or `False` accordingly. The argument is
evaluated first (ordinary, non-held evaluation), so `BooleanQ[2 > 1]` sees the already
reduced `True`. A wrong argument count routes through `builtin_arg_error`.

**Data structures.** Plain pointer comparison against the two interned boolean names;
no allocation beyond the returned result symbol.

**Complexity / limits.** `O(1)`. Attributes: `Protected`. Unlike `TrueQ` — which asks
"does this reduce to `True`?" and so maps everything that is not `True` to `False` —
`BooleanQ` tests the symbol itself, so `BooleanQ[False]` is `True` while `TrueQ[False]`
is `False`. Any non-boolean (a number, an unevaluated symbol, an unreduced
inequality) gives `False`.

**Attributes:** `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`BooleanQ[expr]` gives `True` when `expr` is literally the symbol `True` or `False`,
and `False` for everything else. The argument is evaluated first, so `BooleanQ[2 > 1]`
sees the reduced `True`.

It tests the *symbol*, which is the key difference from `TrueQ`: `BooleanQ[False]` is
`True` (it is a boolean), whereas `TrueQ[False]` is `False` (it is not `True`).
`BooleanQ` is `Protected`.
