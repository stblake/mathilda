# JoinAcross

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`JoinAcross[{a1, ...}, {b1, ...}, spec]`**

Joins two lists of associations, merging each ai with every bj whose join keys agree (inner join). spec is a key, Key\[k\], k1 -\> k2 (keys named differently on the two sides), or a list of these.

**`JoinAcross[{a1, ...}, {b1, ...}, spec, type]`**

type is "Inner", "Left", "Right" or "Outer"; unmatched rows get Missing\["Unmatched"\] for the other side's keys. Option KeyCollisionFunction -\> Left | Right | f resolves a non-join key present on both sides (f\[k\] gives the pair of new keys).

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= JoinAcross[{<|"id" -> 1, "name" -> "Ada"|>, <|"id" -> 2, "name" -> "Bob"|>}, {<|"id" -> 1, "dept" -> "R&D"|>, <|"id" -> 3, "dept" -> "Ops"|>}, Key["id"]]
Out[1]= {<|"id" -> 1, "name" -> "Ada", "dept" -> "R&D"|>}

In[2]:= JoinAcross[{<|"id" -> 1, "name" -> "Ada"|>, <|"id" -> 2, "name" -> "Bob"|>}, {<|"id" -> 1, "dept" -> "R&D"|>, <|"id" -> 3, "dept" -> "Ops"|>}, "id", "Outer"]
Out[2]= {<|"id" -> 1, "name" -> "Ada", "dept" -> "R&D"|>, <|"id" -> 2, "name" -> "Bob", "dept" -> Missing["Unmatched"]|>, <|"id" -> 3, "name" -> Missing["Unmatched"], "dept" -> "Ops"|>}

In[3]:= JoinAcross[{<|"k" -> 1, "x" -> 10|>}, {<|"key" -> 1, "y" -> 20|>}, "k" -> "key"]
Out[3]= {<|"k" -> 1, "x" -> 10, "key" -> 1, "y" -> 20|>}
```

### Options (1)

```mathematica
In[4]:= JoinAcross[{<|"a" -> 1, "v" -> "L"|>}, {<|"a" -> 1, "v" -> "R"|>}, "a", KeyCollisionFunction -> Right]
Out[4]= {<|"a" -> 1, "v" -> "R"|>}
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
