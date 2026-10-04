# KeyMap

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyMap[f, assoc]`**

Applies f to every key, keeping values (keys that collide collapse with last-value-wins).

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= KeyMap[f, <|1 -> 10, 2 -> 20|>]
Out[1]= <|f[1] -> 10, f[2] -> 20|>
```

### Applications (2)

Transform each key, values kept

```mathematica
In[2]:= KeyMap[f, <|a -> 1, b -> 2|>]
Out[2]= <|f[a] -> 1, f[b] -> 2|>
```

Distinct keys may collapse; last value wins

```mathematica
In[3]:= KeyMap[EvenQ, <|1 -> a, 2 -> b, 3 -> c, 4 -> d|>]
Out[3]= <|False -> c, True -> d|>
```

## Implementation notes

**Algorithm.** `builtin_keymap` applies `f` to each key (`apply1(f, key)`), keeps the
value unchanged, and rebuilds the association with `assoc_from_rules`. The
re-canonicalisation matters: `f` may map two distinct keys to the same new key, and
`assoc_from_rules` then applies the usual association rule — a later entry's value
wins for a collided key.

**Data structures.** A fresh `Rule[newkey, value]` array fed to `assoc_from_rules`,
which de-duplicates keys through its own hash index.

**Complexity / limits.** O(n), one `f` evaluation per key, plus the O(n) canonicalise.
Because keys can collide, the result may be shorter than the input. `KeyValueMap`
instead feeds both key and value to `f` and returns a plain `List`.

**Attributes:** `Protected`.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

Only the keys are transformed; the values ride along. If `f` sends two keys to the
same result the association re-canonicalises and the later entry's value survives, so
the output can be shorter than the input (as in the `EvenQ` example, where four keys
collapse to `False` and `True`).
