# Xor

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Xor[e1, e2, ...]`**

The logical exclusive OR of the ei: True when an odd number of the arguments are True.  Folds literal Booleans and cancels duplicate arguments (a Xor a is False); Xor\[\] is False and Xor\[e\] is e.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Xor[True, False]
Out[1]= True

In[2]:= Xor[p, q, p]
Out[2]= q

In[3]:= Xor[True, a]
Out[3]= Not[a]
```

### Applications (5)

An odd number of True arguments gives True

```mathematica
In[4]:= Xor[True, False]
Out[4]= True
```

Three Trues: still odd

```mathematica
In[5]:= Xor[True, True, True]
Out[5]= True
```

Duplicate arguments cancel in pairs

```mathematica
In[6]:= Xor[p, q, p]
Out[6]= q
```

A single True negates the remaining argument

```mathematica
In[7]:= Xor[True, a]
Out[7]= Not[a]
```

Distinct symbolic arguments stay symbolic

```mathematica
In[8]:= Xor[p, q]
Out[8]= Xor[p, q]
```

## Implementation notes

**Algorithm.** `Xor` is `Flat, Orderless, OneIdentity, Protected`, so by the time
`builtin_xor` runs the evaluator has already flattened nested `Xor` calls and
sorted the arguments canonically. The builtin makes one pass over the arguments:
each literal `True` flips a parity bit (and is dropped), each `False` is dropped,
and any argument structurally equal (`expr_eq`) to one already kept cancels
against it (`a` Xor `a` = `False`). What survives is the deduplicated non-literal
core; if an odd number of `True`s were consumed, the result is wrapped in `Not`.
`Xor[]` is `False`, `Xor[e]` collapses to `e` by `OneIdentity`, and when nothing
simplified the builtin returns `NULL` to stay symbolic.

**Data structures.** A single `malloc`'d array of *borrowed* argument pointers
holds the surviving terms; cancelled slots are `NULL`'d and then compacted, and
only the final core is deep-copied into the result. The parity is one `int`.

**Complexity / limits.** The duplicate check is pairwise, so the pass is `O(n²)`
in the argument count with `expr_eq` comparisons — fine for the small boolean
expressions this is used on. Simplification is purely structural: it folds literal
Booleans and exact duplicates but does not reason about implications between
distinct symbolic arguments.

**Attributes:** `Flat`, `OneIdentity`, `Orderless`, `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [Orderless](../../expression-information/Orderless/), [OneIdentity](../../expression-information/OneIdentity/)

- Source: [`src/boolean.c`](https://github.com/stblake/mathilda/blob/main/src/boolean.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_boolean.c`](https://github.com/stblake/mathilda/blob/main/tests/test_boolean.c)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)

## Notes & additional examples

### Notes

`Xor[e1, e2, …]` is `True` when an odd number of the `ei` are `True`. It is
`Flat`, `Orderless` and `OneIdentity`, so nested `Xor` flattens, arguments are
sorted canonically, and `Xor[e]` collapses to `e`.

Evaluation folds the literal Booleans and cancels duplicate arguments (`a` Xor `a`
is `False`): `Xor[]` is `False`, `Xor[True, False]` is `True`, `Xor[True, True]`
is `False`, and an odd count of consumed `True`s negates the surviving core, so
`Xor[True, a]` is `Not[a]`.
