# Select

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Select[list, crit]`**

selects elements e of list for which crit\[e\] yields True, preserving the head of list.

**`Select[list, crit, n]`**

stops after the first n matching elements.

**`Select[crit]`**

is the operator form: Select\[crit\]\[list\] == Select\[list, crit\].

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Map[#^2 &, <|"x" -> 3, "y" -> 4|>]
Out[1]= <|"x" -> 9, "y" -> 16|>

In[2]:= Select[<|"a" -> 1, "b" -> 2, "c" -> 3|>, # > 1 &]
Out[2]= <|"b" -> 2, "c" -> 3|>
```

### Applications (6)

```mathematica
In[3]:= Select[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[3]= {2, 4, 6}

In[4]:= Select[Range[10], # > 5 &]
Out[4]= {6, 7, 8, 9, 10}

In[5]:= Select[{1, 2, 3, 4, 5}, PrimeQ, 2]
Out[5]= {2, 3}

In[6]:= Select[Range[100], PrimeQ[#] && PrimeQ[# + 2] &]
Out[6]= {3, 5, 11, 17, 29, 41, 59, 71}

In[7]:= Select[Range[2, 50], PrimeQ[2^# - 1] &]
Out[7]= {2, 3, 5, 7, 13, 17, 19, 31}

In[8]:= Select[Range[1, 20], GCD[#, 20] == 1 &]
Out[8]= {1, 3, 7, 9, 11, 13, 17, 19}
```

## Implementation notes

`builtin_select` filters the arguments of a compound expression by a predicate.
It iterates the args of `list` (any head, not only `List`), and for each element
builds `crit[elem]` and runs `evaluate()`; the element is kept only when the
result is exactly the symbol `True`. The optional third argument caps the number
of kept elements (`n_max`), stopping the scan early once reached. The surviving
elements are reassembled under the original head via `expr_new_function`. Returns
`NULL` (unevaluated) when the first argument is an atom or when the count
argument is non-integer. Each predicate test allocates a copied call and frees it
plus its evaluated result, so memory is bounded per element.

**Attributes:** `Protected`.

## References

**See also:** [Map](../../data-structures/Map/)

- Harold Abelson and Gerald Jay Sussman, *Structure and Interpretation of Computer Programs*, 2nd ed., §2.2.3 (sequences as conventional interfaces; filtering).
- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
- Tests: [`tests/test_divisors.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisors.c)

## Notes & additional examples

### Notes

`Select[list, crit]` keeps the elements for which `crit[elem]` returns `True`;
any other result (including `False` or an unevaluated predicate) drops the
element. The criterion is usually a predicate symbol such as `EvenQ` or
`PrimeQ`, or a pure function like `# > 5 &`. The optional third argument caps the
number of elements returned, which lets `Select` stop early once enough matches
are found.
