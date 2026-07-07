# GroupBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
GroupBy[list, f]
    Groups the elements of list by the value of f[element],
    giving <|f[x] -> {matching elements}, ...|>.
GroupBy[list, f, g]
    Applies the reducer g to each group, giving
    <|f[x] -> g[{matching elements}], ...|> (split-apply-combine).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= GroupBy[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[1]= <|False -> {1, 3, 5}, True -> {2, 4, 6}|>

In[2]:= GroupBy[Range[10], EvenQ, Total]
Out[2]= <|False -> 25, True -> 30|>

In[3]:= GroupBy[{{"x", 1}, {"y", 2}, {"x", 3}}, First -> Last, Total]
Out[3]= <|"x" -> 4, "y" -> 2|>

In[4]:= GroupBy[<|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4|>, EvenQ]
Out[4]= <|False -> <|"a" -> 1, "c" -> 3|>, True -> <|"b" -> 2, "d" -> 4|>|>

In[5]:= GroupBy[<|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4|>, EvenQ, Total]
Out[5]= <|False -> 4, True -> 6|>
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
