# Discard

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Discard[expr, crit]`**

Drops the elements e of expr for which crit\[e\] is True (the complement of Select). On an association crit tests the values.

**`Discard[expr, crit, n]`**

Drops only the first n such elements (n may be Infinity).

## Examples (7)

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

### Applications (4)

The complement of Select

```mathematica
In[4]:= Discard[{1, 2, 3, 4, 5}, EvenQ]
Out[4]= {1, 3, 5}
```

```mathematica
In[5]:= Discard[Range[10], # > 5 &]
Out[5]= {1, 2, 3, 4, 5}
```

Drop at most the first match

```mathematica
In[6]:= Discard[{1, 2, 3, 4, 5}, OddQ, 1]
Out[6]= {2, 3, 4, 5}
```

Over an association, tests the values

```mathematica
In[7]:= Discard[<|a -> 1, b -> 2, c -> 3|>, OddQ]
Out[7]= <|b -> 2|>
```

## Implementation notes

**Algorithm.** `builtin_discard` is the complement of `Select`: it keeps the
elements of a non-atomic expression for which the criterion does *not* give
`True`. `Discard[expr, crit]` tests every element; `Discard[expr, crit, n]`
discards at most the first `n` matches and keeps the rest regardless. The result
reuses the original head, so it works on any non-atomic expression, and on an
association it tests the values while keeping the keys.

**Data structures.** A single pass builds a `keep` array of element copies; the
criterion is applied with `eval_call1` (evaluate `crit[elem]`) and the verdict
read with `is_true_sym`. A visible `NDArray` is de-listed first via
`ops_delist_visible`.

**Complexity / limits.** `O(n)` evaluations of the criterion. A non-atomic first
argument is required (`Discard::normal`), and the optional limit must be a
non-negative integer or `Infinity` (`Discard::innf`); both diagnostics route
through the message funnel so `Quiet[]`/`Check[]` behave. `crit` is arbitrary, so
there is no packed fast path.

**Attributes:** `Protected`.

## References

**See also:** [Select](../../data-structures/Select/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)

## Notes & additional examples

### Notes

`Discard[expr, crit]` keeps the elements for which `crit` does *not* return
`True` — exactly the complement of `Select[expr, crit]`. The optional third
argument `Discard[expr, crit, n]` discards at most the first `n` matching
elements and keeps everything after. It works on any non-atomic expression and,
over an association, tests the values while preserving keys.
