# Lookup

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Lookup[assoc, key]`**

Gives the value for key, or Missing\["KeyAbsent", key\].

**`Lookup[assoc, key, default]`**

Uses default when key is absent.

**`Lookup[assoc, {k1, k2, ...}]`**

Looks up several keys at once (O(n+m)).

**`Lookup[{assoc1, assoc2, ...}, key, ...]`**

Threads over a list of associations and/or lists of rules.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Lookup[<|"a" -> 1, "b" -> 2|>, "b"]
Out[1]= 2

In[2]:= Lookup[<|"a" -> 1|>, "z", 0]
Out[2]= 0

In[3]:= Lookup[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3|>}, "a", 0]
Out[3]= {1, 3}

In[4]:= Lookup[<|"a" -> 1|>, "a", Print["never printed"]; 0]
Out[4]= 1

In[5]:= Lookup["a"][<|"a" -> 7|>]
Out[5]= 7
```

## Implementation notes

- `HoldAll`, as in Mathematica: the association and key are evaluated, the
  default only when a key is actually absent (once per absent key), so a
  default with side effects runs only when needed.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Keys](../../data-structures/Keys/), [Values](../../data-structures/Values/), [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
