# KeyMemberQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyMemberQ[assoc, patt]`**

Gives True if some key of assoc matches the pattern patt (a literal key is an O(1) probe), else False.

## Examples (8)

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

### Applications (3)

A is a key

```mathematica
In[6]:= KeyMemberQ[<|a -> 1, b -> 2|>, a]
Out[6]= True
```

C is not

```mathematica
In[7]:= KeyMemberQ[<|a -> 1, b -> 2|>, c]
Out[7]= False
```

The second argument is a pattern

```mathematica
In[8]:= KeyMemberQ[<|1 -> x, 2 -> y|>, _Integer]
Out[8]= True
```

## Implementation notes

**Algorithm.** `builtin_keymemberq` asks whether any key of the association matches
its second argument, which is read as a *pattern*. `assoc_some_key_matches` first
calls `key_contains_pattern`: a pattern-free query is answered by a single
`assoc_lookup_value` probe (O(1) index lookup); a query carrying a pattern
(`Blank`, `Pattern`, `Alternatives`, `PatternTest`, …) is matched structurally
against each key, short-circuiting on the first hit. An association or a bare list
of rules is accepted.

**Data structures.** The association's persistent open-addressing hash index
(`KeyIndex`) for the literal path; a transient `MatchEnv` per key for the pattern
path.

**Complexity / limits.** O(1) amortised for a literal key, O(n) for a pattern.
Unlike `KeyExistsQ` (literal argument), `KeyMemberQ[a, _]` is `True` for any
non-empty `a`.

- A pattern-free key keeps the O(1) index probe; a key containing a pattern
  construct (`_`, `x_h`, `Alternatives`, `Except`, `PatternTest`, `Condition`,
  ...) is matched against every key.
- `KeyExistsQ` never treats its key as a pattern: `KeyExistsQ[a, _]` looks for
  the literal key `_`, as in Mathematica.

**Attributes:** `Protected`.

## References

**See also:** [KeyExistsQ](../../data-structures/KeyExistsQ/), [KeyFreeQ](../../data-structures/KeyFreeQ/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)

## Notes & additional examples

### Notes

The second argument is a **pattern**: `KeyMemberQ[a, _]` is `True` for any non-empty
association, whereas `KeyExistsQ[a, _]` searches for the literal key `_`. A
pattern-free key takes the O(1) hash-index probe; a pattern falls back to a match
over each key. `KeyFreeQ` is the complement.
