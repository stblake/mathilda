# GroupBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GroupBy[list, f]`**

Groups the elements of list by the value of f\[element\], giving \<|f\[x\] -\> {matching elements}, ...|\>.

**`GroupBy[list, f, g]`**

Applies the reducer g to each group, giving \<|f\[x\] -\> g\[{matching elements}\], ...|\> (split-apply-combine).

## Examples (11)

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

### Applications (4)

```mathematica
In[8]:= GroupBy[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[8]= <|False -> {1, 3, 5}, True -> {2, 4, 6}|>

In[9]:= GroupBy[Range[6], Mod[#, 3] &]
Out[9]= <|1 -> {1, 4}, 2 -> {2, 5}, 0 -> {3, 6}|>
```

A reducer summarises each group

```mathematica
In[10]:= GroupBy[{1, 2, 3, 4, 5, 6}, EvenQ, Total]
Out[10]= <|False -> 9, True -> 12|>
```

```mathematica
In[11]:= GroupBy[{"apple", "pear", "plum", "fig"}, StringLength]
Out[11]= <|5 -> {"apple"}, 4 -> {"pear", "plum"}, 3 -> {"fig"}|>
```

## Implementation notes

**Algorithm.** `builtin_groupby` returns `<|f[x] -> {x, ...}|>`, preserving
first-key order. It evaluates the key function once per element, looks the key up
in a `KeyIndex`, and appends the element (a copy) to that key's growing group
buffer. A third argument reduces each group (`GroupBy[list, f, red]` →
`<|k -> red[group]|>`); a `keyfn -> valfn` second argument groups by `keyfn[x]`
but collects `valfn[x]`; and a `List` of key functions grooves multi-level
grouping into nested associations, the reducer applying only at the innermost
level. Over an `Association` the entries are grouped by `f[value]` with keys
preserved.

**Data structures.** A `KeyIndex` open-addressing hash set over the owned group
keys, parallel with per-group dynamic `Expr**` buffers that double on demand; the
groups become `List`s (or sub-associations) under the group keys.

**Complexity / limits.** `O(n)` evaluations of the key function plus `O(n)`
hashing and copying; multi-level grouping recurses through the evaluator (each
level a `GroupBy` call). Returns `NULL` unless the first argument is a `List` or
association. `f` is an arbitrary function, so there is no packed fast path.

**Attributes:** `Protected`.

## References

**See also:** [Total](../../arithmetic/Total/)

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

`GroupBy[list, f]` returns `<|f[x] -> {elements}, ...|>`, grouping the elements by
the value of the key function and preserving first-key order. The three-argument
form `GroupBy[list, f, red]` applies a reducer to each group — `GroupBy[data,
key, Total]` is a group-and-sum in one step. A `keyfn -> valfn` second argument
groups by one function but collects another, and a list of key functions groups
into nested associations. It is the keyed counterpart of `GatherBy`; over an
association the entries are grouped by `f[value]` with keys preserved.
