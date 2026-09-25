# GroupBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GroupBy[list, f]`**

Groups the elements of list by the value of f\[element\], giving \<|f\[x\] -\> {matching elements}, ...|\>.

**`GroupBy[list, f, g]`**

Applies the reducer g to each group, giving \<|f\[x\] -\> g\[{matching elements}\], ...|\> (split-apply-combine).

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= GroupBy[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[1]= <|False -> {1, 3, 5}, True -> {2, 4, 6}|>

In[2]:= GroupBy[Range[10], EvenQ, Total]
Out[2]= <|False -> 25, True -> 30|>

In[3]:= GroupBy[{1, 2, 3, 4, 5, 6}, {EvenQ, # > 3 &}]
Out[3]= <|False -> <|False -> {1, 3}, True -> {5}|>, True -> <|False -> {2}, True -> {4, 6}|>|>

In[4]:= GroupBy[{1, 2, 3, 4, 5, 6}, {EvenQ, # > 3 &}, Total]
Out[4]= <|False -> <|False -> 4, True -> 5|>, True -> <|False -> 2, True -> 10|>|>

In[5]:= GroupBy[<|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4|>, EvenQ]
Out[5]= <|False -> <|"a" -> 1, "c" -> 3|>, True -> <|"b" -> 2, "d" -> 4|>|>

In[6]:= GroupBy[<|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4|>, EvenQ, Total]
Out[6]= <|False -> 4, True -> 6|>
```

### Options (1)

```mathematica
In[7]:= GroupBy[{{"x", 1}, {"y", 2}, {"x", 3}}, First -> Last, Total]
Out[7]= <|"x" -> 4, "y" -> 2|>
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Total](../../arithmetic/Total/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
