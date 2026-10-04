# KeyFreeQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyFreeQ[assoc, patt]`**

Gives True if no key of assoc matches the pattern patt (the complement of KeyMemberQ), else False.

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

True: the key c is absent

```mathematica
In[6]:= KeyFreeQ[<|a -> 1, b -> 2|>, c]
Out[6]= True
```

False: a is present

```mathematica
In[7]:= KeyFreeQ[<|a -> 1, b -> 2|>, a]
Out[7]= False
```

A pattern: some key matches, so not free

```mathematica
In[8]:= KeyFreeQ[<|1 -> x, 2 -> y|>, _Integer]
Out[8]= False
```

## Implementation notes

**Algorithm.** `builtin_keyfreeq` is the boolean complement of `KeyMemberQ`: both
call `assoc_some_key_matches`, and `KeyFreeQ` negates its verdict. The second
argument is treated as a *pattern*, as in Mathematica. When it is pattern-free
(`key_contains_pattern` finds no `Blank`/`Pattern`/`Alternatives`/… construct) the
test is a single `assoc_lookup_value` probe — the O(1) index lookup. When it does
carry a pattern, the matcher is run against each key in turn (`match` into a fresh
`MatchEnv`) and the first match wins. Either an association or a bare list of rules
is accepted.

**Data structures.** For the literal-key path, the association's persistent
open-addressing hash index (`KeyIndex`, keys compared by `expr_eq` and hashed by
`expr_hash`); for the pattern path, a throwaway `MatchEnv` per key.

**Complexity / limits.** O(1) amortised for a literal key, O(n) for a pattern (one
structural match per entry). Contrast `KeyExistsQ`, which looks its argument up
literally even when it is a `_`.

- A pattern-free key keeps the O(1) index probe; a key containing a pattern
  construct (`_`, `x_h`, `Alternatives`, `Except`, `PatternTest`, `Condition`,
  ...) is matched against every key.
- `KeyExistsQ` never treats its key as a pattern: `KeyExistsQ[a, _]` looks for
  the literal key `_`, as in Mathematica.

**Attributes:** `Protected`.

## References

**See also:** [KeyExistsQ](../../data-structures/KeyExistsQ/), [KeyMemberQ](../../data-structures/KeyMemberQ/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)

## Notes & additional examples

### Notes

`KeyFreeQ` is the exact complement of `KeyMemberQ`. Its second argument is a
**pattern**, so `KeyFreeQ[a, _]` is `False` for any non-empty association. This is
the one place it differs from `KeyExistsQ`, which looks for the literal key `_`.
Both an association and a bare list of rules are accepted.
