# Association

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Association[key1 -> val1, key2 -> val2, ...]  (also written <|...|>)`**

Represents an association mapping keys to values with unique, insertion-ordered keys (last value wins on duplicates). Arguments may be rules, lists of rules, or other associations.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= <|"a" -> 1, "b" -> 2|>
Out[1]= <|"a" -> 1, "b" -> 2|>

In[2]:= <|"a" -> 1, "b" -> 2, "a" -> 99|>
Out[2]= <|"a" -> 99, "b" -> 2|>

In[3]:= <|"a" -> 10, "b" -> 20|>[["b"]]
Out[3]= 20

In[4]:= <|"a" -> 10, "b" -> 20|>["a"]
Out[4]= 10

In[5]:= <|"a" -> <|"b" -> 5|>|>["a", "b"]
Out[5]= 5

In[6]:= Association[{{"a" -> 1}, {"b" -> 2, {"c" -> 3}}}]
Out[6]= <|"a" -> 1, "b" -> 2, "c" -> 3|>
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_core_algebra.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core_algebra.c)
