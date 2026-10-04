# KeyExistsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyExistsQ[assoc, key]`**

Gives True if key is present in assoc, else False.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= KeyExistsQ[<|"a" -> 1|>, "a"]
Out[1]= True

In[2]:= KeyFreeQ[<|"a" -> 1|>, "b"]
Out[2]= True

In[3]:= KeyMemberQ[<|1 -> "x", "b" -> "y"|>, _Integer]
Out[3]= True

In[4]:= KeyFreeQ[<|"a" -> 1|>, _Integer]
Out[4]= True

In[5]:= KeyExistsQ[<|"a" -> 1|>, _]
Out[5]= False
```

### Applications (2)

```mathematica
In[6]:= KeyExistsQ[<|a -> 1, b -> 2|>, a]
Out[6]= True

In[7]:= KeyExistsQ[<|a -> 1, b -> 2|>, c]
Out[7]= False
```

## Implementation notes

**Algorithm.** `builtin_keyexistsq` returns `True` iff the literal key is present
in the association (or bare list of rules). It is a single `assoc_lookup_value`
probe: a non-`NULL` value means the key exists. Unlike `KeyMemberQ`/`KeyFreeQ`,
the second argument is treated as a *literal* key, not a pattern — so
`KeyExistsQ[a, _]` looks for the key spelled `_`, not "any key".

**Data structures.** `assoc_lookup_value` uses the association's cached
`AssocIndex` when present (built lazily on the first single-key read), giving an
`O(1)` open-addressing hash probe; otherwise it falls back to an `O(n)`
`assoc_scan`. Keys are compared with `expr_eq`.

**Complexity / limits.** `O(1)` amortised once the key index exists, else `O(n)`.
Returns `NULL` (unevaluated) unless the first argument is an association or a list
of rules; it is a `*Q` predicate and otherwise always answers a Boolean.

- A pattern-free key keeps the O(1) index probe; a key containing a pattern
  construct (`_`, `x_h`, `Alternatives`, `Except`, `PatternTest`, `Condition`,
  ...) is matched against every key.
- `KeyExistsQ` never treats its key as a pattern: `KeyExistsQ[a, _]` looks for
  the literal key `_`, as in Mathematica.

**Attributes:** `Protected`.

## References

**See also:** [KeyMemberQ](../../data-structures/KeyMemberQ/), [KeyFreeQ](../../data-structures/KeyFreeQ/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
- Tests: [`tests/test_parallelmixedtower.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parallelmixedtower.c)

## Notes & additional examples

### Notes

`KeyExistsQ[assoc, key]` returns `True` if the association has the given key, else
`False`. The key is taken *literally*, not as a pattern — this is what
distinguishes `KeyExistsQ` from `KeyMemberQ`/`KeyFreeQ`, which match their second
argument as a pattern (`KeyMemberQ[a, _]` is `True` for any non-empty `a`, while
`KeyExistsQ[a, _]` looks for the literal key `_`). The test is a single `O(1)`
amortised probe through the association's key index, and it also accepts a bare
list of rules.
