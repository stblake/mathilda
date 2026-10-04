# Delete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Delete[expr, n] deletes the element at position n in expr.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Delete[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {Key["b"]}]
Out[1]= <|"a" -> 1, "c" -> 3|>

In[2]:= Delete[<|"a" -> <|"x" -> 5, "y" -> 6|>|>, {Key["a"], Key["x"]}]
Out[2]= <|"a" -> <|"y" -> 6|>|>
```

### Applications (4)

```mathematica
In[3]:= Delete[{a, b, c, d}, 2]
Out[3]= {a, c, d}
```

Negative indices count from the end

```mathematica
In[4]:= Delete[{a, b, c, d}, -1]
Out[4]= {a, b, c}
```

A list of positions deletes several at once

```mathematica
In[5]:= Delete[{a, b, c, d}, {{1}, {3}}]
Out[5]= {b, d}
```

Works on any head

```mathematica
In[6]:= Delete[f[a, b, c], 2]
Out[6]= f[a, c]
```

## Implementation notes

`builtin_delete` (in `src/part.c`) drives the recursive helper `delete_path`, which walks an integer position (or position path) into the expression tree and removes the targeted element by rebuilding the enclosing function with all arguments except that index. Negative indices count from the end, position `0` targets the head (replaced by a `Sequence[...]` of the remaining parts), and out-of-range indices leave the structure unchanged.

**Attributes:** none registered.

## References

**See also:** [KeyDrop](../../data-structures/KeyDrop/)

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_part.c`](https://github.com/stblake/mathilda/blob/main/tests/test_part.c)

## Notes & additional examples

### Notes

`Delete[expr, n]` removes the element at position `n`, rebuilding the enclosing
expression without it; negative indices count from the end, and a position path
`{i, j, ...}` reaches a nested element. A list of positions
`{{p1}, {p2}, ...}` deletes each of them in one call. `Delete` works on any head,
not just `List`. Unlike `Drop`, which removes a count or a strided range, `Delete`
targets positions explicitly; out-of-range positions leave the structure
unchanged.
