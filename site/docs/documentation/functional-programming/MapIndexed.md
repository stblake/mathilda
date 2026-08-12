# MapIndexed

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MapIndexed[f, expr]`**

Applies f to the elements of expr, giving the part specification of each element as a second argument to f: {f\[e1, {1}\], f\[e2, {2}\], ...}.

**`MapIndexed[f, expr, levelspec]`**

Applies f to all parts of expr on the levels specified by levelspec: n            levels 1 through n Infinity     levels 1 through Infinity {n}          level n only {n1, n2}     levels n1 through n2 The default is {1}. A positive level n consists of all parts specified by n indices; a negative level -n consists of all parts with depth n, so level -1 is the atoms. Level 0 is the whole expression, whose position is {}. The position handed to f is the one Part and Extract take, so Extract\[expr, #2\] is #1. Over an association the position of a value is {Key\[k\]}, and keys are preserved. With Heads -\> True the function is applied to heads as well, a head having index 0 in its position. MapIndexed always effectively constructs a complete new expression and then evaluates it.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= MapIndexed[f, {10, 20, 30}]
Out[1]= {f[10, {1}], f[20, {2}], f[30, {3}]}

In[2]:= MapIndexed[First[#2] + f[#1] &, {a, b, c, d}]
Out[2]= {1 + f[a], 2 + f[b], 3 + f[c], 4 + f[d]}

In[3]:= MapIndexed[f, {{a, b}, {c, d, e}}, {2}]
Out[3]= {{f[a, {1, 1}], f[b, {1, 2}]}, {f[c, {2, 1}], f[d, {2, 2}], f[e, {2, 3}]}}

In[4]:= MapIndexed[f, {{a, b}, {c, d, {e}}}, {-1}]
Out[4]= {{f[a, {1, 1}], f[b, {1, 2}]}, {f[c, {2, 1}], f[d, {2, 2}], {f[e, {2, 3, 1}]}}}

In[5]:= MapIndexed[f, h0[h1[h2[h3[h4[a]]]]], {2, -3}]
Out[5]= h0[h1[f[h2[f[h3[h4[a]], {1, 1, 1}]], {1, 1}]]]

In[6]:= MapIndexed[f, {a, b}, {0, 1}]
Out[6]= f[{f[a, {1}], f[b, {2}]}, {}]

In[7]:= MapIndexed[f, <|"a" -> 10, "b" -> 20|>]
Out[7]= <|"a" -> f[10, {Key["a"]}], "b" -> f[20, {Key["b"]}]|>
```

### Options (1)

```mathematica
In[8]:= MapIndexed[f, p[x][a, b, c], Infinity, Heads -> True]
Out[8]= f[f[p, {0, 0}][f[x, {0, 1}]], {0}][f[a, {1}], f[b, {2}], f[c, {3}]]
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[Part](../../structural-manipulation/Part/), [Extract](../../structural-manipulation/Extract/), [Rule](../../assignment-and-rules/Rule/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/), [Association](../../data-structures/Association/), [List](../../other-advanced/List/)

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_map_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_map_ndarray.c)
- Tests: [`tests/test_mapindexed.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mapindexed.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
