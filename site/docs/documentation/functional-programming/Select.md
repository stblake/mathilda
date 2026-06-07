# Select

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Select[list, crit]
    selects elements e of list for which crit[e] yields True, preserving
    the head of list.
Select[list, crit, n]
    stops after the first n matching elements.
Select[crit]
    is the operator form: Select[crit][list] == Select[list, crit].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Select[{1, 2, 4, 7, 6, 2}, EvenQ]
Out[1]= {2, 4, 6, 2}

In[2]:= Select[{1, 2, 4, 7, 6, 2}, # > 2 &, 1]
Out[2]= {4}
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

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Harold Abelson and Gerald Jay Sussman, *Structure and Interpretation of Computer Programs*, 2nd ed., §2.2.3 (sequences as conventional interfaces; filtering).
- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Select[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[1]= {2, 4, 6}
```

```mathematica
In[1]:= Select[Range[10], # > 5 &]
Out[1]= {6, 7, 8, 9, 10}
```

```mathematica
In[1]:= Select[{1, 2, 3, 4, 5}, PrimeQ, 2]
Out[1]= {2, 3}
```

### Notes

`Select[list, crit]` keeps the elements for which `crit[elem]` returns `True`;
any other result (including `False` or an unevaluated predicate) drops the
element. The criterion is usually a predicate symbol such as `EvenQ` or
`PrimeQ`, or a pure function like `# > 5 &`. The optional third argument caps the
number of elements returned, which lets `Select` stop early once enough matches
are found.
