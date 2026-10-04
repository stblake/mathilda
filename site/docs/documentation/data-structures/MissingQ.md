# MissingQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MissingQ[expr]`**

Gives True if expr has head Missing (Missing\[\], Missing\["reason"\], Missing\["KeyAbsent", k\], ...), and False otherwise.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= MissingQ[Missing["KeyAbsent", "z"]]
Out[1]= True

In[2]:= MissingQ[<|"a" -> 1|>["z"]]
Out[2]= True

In[3]:= Select[{1, Missing[], 3}, Not @* MissingQ]
Out[3]= {1, 3}
```

### Applications (3)

```mathematica
In[4]:= MissingQ[Missing[]]
Out[4]= True
```

The absent-key result is Missing

```mathematica
In[5]:= MissingQ[Lookup[<|a -> 1|>, b]]
Out[5]= True
```

```mathematica
In[6]:= MissingQ[5]
Out[6]= False
```

## Implementation notes

**Algorithm.** `builtin_missingq` returns `True` exactly when its single argument has
head `Missing` — `head_is(arg, SYM_Missing)` — for any arity (`Missing[]`,
`Missing["reason"]`, `Missing["reason", data]`). A call with other than one argument
emits `MissingQ::argx` and is left unevaluated.

**Data structures.** None beyond the head comparison against the interned
`SYM_Missing`.

**Complexity / limits.** O(1). It is purely a head test, so it says nothing about
*why* the value is missing; pair it with `DeleteMissing` to remove such entries, or
with `Lookup`'s default to replace them.

**Attributes:** `Protected`.

## References

**See also:** [Missing](../../data-structures/Missing/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)

## Notes & additional examples

### Notes

`MissingQ` is a plain head test: `True` for any `Missing[…]` expression and `False`
otherwise. It is the standard guard after a `Lookup` or key access that may not find
its key.
