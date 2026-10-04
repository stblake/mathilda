# AssociateTo

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AssociateTo[s, key -> val]  |  AssociateTo[s, {rules}]`**

Adds or updates key-value pairs in the association held by symbol s, modifying s in place.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= asc = <|"a" -> 1|>; AssociateTo[asc, "b" -> 2]; asc
Out[1]= <|"a" -> 1, "b" -> 2|>
```

### Applications (4)

```mathematica
In[2]:= a = <|x -> 1, y -> 2|>;
```

A new key is appended at the end

```mathematica
In[3]:= AssociateTo[a, z -> 3]
Out[3]= <|x -> 1, y -> 2, z -> 3|>
```

An existing key is updated in place

```mathematica
In[4]:= AssociateTo[a, x -> 9]
Out[4]= <|x -> 9, y -> 2, z -> 3|>
```

```mathematica
In[5]:= a
Out[5]= <|x -> 9, y -> 2, z -> 3|>
```

## Implementation notes

**Algorithm.** `builtin_associate_to` is the in-place update for associations,
mirroring `AppendTo`. It is `HoldFirst`, so the first argument is the unevaluated
symbol; the builtin evaluates it to read the symbol's current association,
gathers the existing entries together with the new rule(s) (a single
`key -> val` or a `List` of rules), re-canonicalises the whole set with
`assoc_from_rules`, assigns the updated association back to the symbol via
`symtab_add_own_value`, and returns it.

**Data structures.** A flat `Expr*` array of the base entries plus the new ones
is passed to `assoc_from_rules`, whose transient `KeyIndex` hash set applies the
*first position, last value* de-duplication — so associating an existing key
overwrites its value in place while a new key is appended at the end.

**Complexity / limits.** `O(n + k)` for `n` existing entries and `k` new rules.
Returns `NULL` (unevaluated) if the first argument is not a symbol, if its
current value is not an association, or if any supplied entry is not a 2-arg
rule.

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [HoldFirst](../../other-advanced/HoldFirst/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)

## Notes & additional examples

### Notes

`AssociateTo[s, key -> val]` is the in-place analogue of `AppendTo` for
associations: it reads the current association stored in `s`, inserts or updates
the entry, assigns the result back to `s`, and returns it. The second argument
may also be a list of rules. Because the key set is de-duplicated with the usual
*first position, last value* rule, associating an existing key overwrites its
value without moving it, while a genuinely new key is added at the end. `s` must
already hold an association.
