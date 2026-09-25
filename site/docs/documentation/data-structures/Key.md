# Key

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Key[k]`**

Represents the key k of an association. Key\[k\]\[assoc\] gives the value at k (Missing\["KeyAbsent", k\] if absent); assoc\[\[Key\[k\]\]\], Lookup and the key-spec arguments of GroupBy, SortBy and JoinAcross accept it, and Key\[k\] names a key literally even when k is an integer or a list.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Key["b"][<|"a" -> 1, "b" -> 2|>]
Out[1]= 2

In[2]:= Lookup[<|1 -> "one", {1, 2} -> "pair"|>, Key[{1, 2}]]
Out[2]= "pair"

In[3]:= Key["z"][<|"a" -> 1|>]
Out[3]= Missing["KeyAbsent", "z"]
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [GroupBy](../../data-structures/GroupBy/), [SortBy](../../data-structures/SortBy/), [JoinAcross](../../data-structures/JoinAcross/), [MapAt](../../data-structures/MapAt/), [Extract](../../structural-manipulation/Extract/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
