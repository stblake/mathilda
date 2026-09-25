# Discard

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Discard[expr, crit]`**

Drops the elements e of expr for which crit\[e\] is True (the complement of Select). On an association crit tests the values.

**`Discard[expr, crit, n]`**

Drops only the first n such elements (n may be Infinity).

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Discard[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[1]= {1, 3, 5}

In[2]:= Discard[<|"a" -> 1, "b" -> 2, "c" -> 4|>, EvenQ]
Out[2]= <|"a" -> 1|>

In[3]:= Discard[{1, 2, 3, 4, 5, 6}, EvenQ, 2]
Out[3]= {1, 3, 5, 6}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Select](../../data-structures/Select/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
