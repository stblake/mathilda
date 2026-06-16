# Outer

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Outer[f,list1,list2,...]
    gives the generalized outer product of the listi, forming all possible combinations of the lowest-level elements in each of them, and feeding them as arguments to f.
Outer[f,list1,list2,...,n]
    treats as separate elements only sublists at level n in the listi.
Outer[f,list1,list2,...,n1,n2,...]
    treats as separate elements only sublists at level ni in the corresponding listi.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Outer[f, {a, b}, {x, y, z}]
Out[1]= {{f[a, x], f[a, y], f[a, z]}, {f[b, x], f[b, y], f[b, z]}}

In[2]:= Outer[Times, {1, 2, 3, 4}, {a, b, c}]
Out[2]= {{a, b, c}, {2 a, 2 b, 2 c}, {3 a, 3 b, 3 c}, {4 a, 4 b, 4 c}}

In[3]:= Outer[g, f[a, b], f[x, y, z]]
Out[3]= f[f[g[a, x], g[a, y], g[a, z]], f[g[b, x], g[b, y], g[b, z]]]
```

## Implementation notes

**Algorithm.** `builtin_outer` computes the generalised outer product `Outer[f, t1, t2, ..., {n1, ...}]`. It first counts trailing Integer / `Infinity` arguments as per-tensor depth limits (default `INT64_MAX`, i.e. descend to the leaves), leaving the remaining arguments after `f` as the input tensors. The recursive worker `outer_rec` walks each tensor down to its target depth, collecting one atom per tensor into `current_atoms`, and at the deepest level emits `f[a1, a2, ...]`; the assembled tree is then `evaluate`-d once. The result head for the assembled levels is taken from the first function-typed tensor.

**Data structures / limits.** Nested `Expr*` `List`s; per-tensor depths in an `int64_t[]`. With no tensors, `f[]` is returned evaluated. This is the generic functional-programming `Outer`, not a linear-algebra-specific kernel; `KroneckerProduct` and matrix outer products are built on it.

- `Protected`.
- Applying `Outer` to two tensors of ranks $r$ and $s$ gives a tensor of rank $r+s$.
- The heads of all `listi` must be the same, but need not necessarily be `List`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

The generalized outer product — every pairing of elements fed to `f`:

```mathematica
In[1]:= Outer[Times, {1, 2, 3}, {a, b, c}]
Out[1]= {{a, b, c}, {2 a, 2 b, 2 c}, {3 a, 3 b, 3 c}}
```

With a symbolic `f` you get the full table of combinations:

```mathematica
In[1]:= Outer[f, {1, 2}, {x, y}]
Out[1]= {{f[1, x], f[1, y]}, {f[2, x], f[2, y]}}
```

`Outer[Times, v, v]` builds the rank-one matrix `v v^T`. Taking the basis
`{1, x, x^2}` produces the symmetric Hankel-style table of monomial products —
the kind of structure that underlies Vandermonde and Gram constructions:

```mathematica
In[1]:= Outer[Times, {1, x, x^2}, {1, x, x^2}]
Out[1]= {{1, x, x^2}, {x, x^2, x^3}, {x^2, x^3, x^4}}
```

Every such outer product of two vectors is rank one, so its determinant must
vanish — a fact the symbolic linear algebra confirms exactly:

```mathematica
In[1]:= Det[Outer[Times, {a, b, c}, {1, 1, 1}]]
Out[1]= 0
```

Outer products of three or more lists nest one level deeper per list:

```mathematica
In[1]:= Outer[List, {1, 2}, {a, b}, {X, Y}]
Out[1]= {{{{1, a, X}, {1, a, Y}}, {{1, b, X}, {1, b, Y}}}, {{{2, a, X}, {2, a, Y}}, {{2, b, X}, {2, b, Y}}}}
```

### Notes

`Outer[f, l1, l2, ...]` forms every combination of lowest-level elements, one
from each list, and applies `f` to it. The result has nesting depth equal to the
combined depth of the inputs, so `Outer` of two vectors is a matrix, of three is
a rank-3 array, and so on. With `f = Times` it is the tensor (outer) product;
with a symbolic head it is a complete combination table. The optional level
arguments `Outer[f, l1, l2, ..., n]` (or per-list `n1, n2, ...`) control which
sublists are treated as the separate elements to combine.
