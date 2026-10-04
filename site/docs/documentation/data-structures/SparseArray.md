# SparseArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SparseArray[{pos1 -> v1, ...}], SparseArray[rules, dims], SparseArray[rules, dims, default]`**

specify an array whose listed positions hold vi and whose other entries are the default (0). Positions are Integer lists (an Integer for a vector); rules may be {p1, p2} -\> {v1, v2}, patterns ({i\_, i\_} -\> 1, needing dims) or Band\[start\] -\> v. SparseArray\[list\] and the internal form SparseArray\[Automatic, dims, default, {1, {rowptr, colidx}, vals}\] are accepted.

<details>
<summary>Notes</summary>

Mathilda has no sparse storage: SparseArray stays inert, and Normal converts it to the dense List.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Normal[SparseArray[{1 -> 1, 3 -> 2}]]
Out[1]= {1, 0, 2}

In[2]:= Normal[SparseArray[{{i_, i_} -> 1}, {3, 3}]]
Out[2]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}

In[3]:= Normal[SparseArray[{{1, 1} -> a, {2, 3} -> b}, {2, 3}, x]]
Out[3]= {{a, x, x}, {x, x, b}}

In[4]:= Normal[SparseArray[Band[{1, 2}] -> {p, q}, {3, 3}]]
Out[4]= {{0, p, 0}, {0, 0, q}, {0, 0, 0}}
```

### Applications (7)

Stays inert: there is no sparse storage

```mathematica
In[5]:= SparseArray[{1 -> a, 3 -> c, 5 -> e}, 5]
Out[5]= SparseArray[{1 -> a, 3 -> c, 5 -> e}, 5]
```

Confirming it does not evaluate

```mathematica
In[6]:= Head[SparseArray[{1 -> 5}, 3]]
Out[6]= SparseArray
```

Normal expands to the dense list

```mathematica
In[7]:= Normal[SparseArray[{1 -> a, 3 -> c, 5 -> e}, 5]]
Out[7]= {a, 0, c, 0, e}
```

A 2x2 identity

```mathematica
In[8]:= Normal[SparseArray[{{1, 1} -> 1, {2, 2} -> 1}, {2, 2}]]
Out[8]= {{1, 0}, {0, 1}}
```

Pattern rule: the 3x3 identity

```mathematica
In[9]:= Normal[SparseArray[{i_, i_} -> 1, {3, 3}]]
Out[9]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
```

A non-zero default

```mathematica
In[10]:= Normal[SparseArray[{2 -> 7}, 4, -1]]
Out[10]= {-1, 7, -1, -1}
```

Band fills a diagonal

```mathematica
In[11]:= Normal[SparseArray[Band[{1, 1}] -> 1, {3, 3}]]
Out[11]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
```

## Implementation notes

**Algorithm.** `SparseArray` is an **inert head** — Mathilda has no sparse storage.
There is no `builtin_sparsearray`: a `SparseArray[rules, dims, default, …]`
specification simply persists unevaluated, and `Head[SparseArray[{1 -> 5}, 3]]` is
`SparseArray`. All the real work lives in `sparse_array_to_dense`, which is reached
*only* through `Normal[SparseArray[…]]` (dispatched in `builtin_normal`,
`src/calculus/series.c`). It allocates one flat row-major buffer of `prod(dims)`
slots, fills it rule by rule (first writer wins), tops the untouched slots up with the
default (`0` unless given), and folds the buffer into nested `List`s. Accepted spec
forms: explicit position rules `{pos -> v}`, `{p1, p2} -> {v1, v2}`, pattern rules
such as `{i_, i_} -> 1` (which visit every slot and so need `dims`), `Band[start] ->
v`, a dense `List`, and the internal CSR form
`SparseArray[Automatic, dims, default, {1, {rowptr, colidx}, vals}]`.

**Data structures.** A transient dense `Dense{rank, dims[32], total, slots}` buffer of
`Expr*` (NULL = not yet written); no compressed/CSR representation is *kept* — the
object remains its symbolic specification. Positions are 1-based with negative
indices resolved from the end (`dense_offset`).

**Complexity / limits.** The `Normal` conversion is O(prod(dims)) in both time and
space, which is the honest cost: because there is no sparse backing store, arithmetic
and linear algebra are **not** accelerated on a `SparseArray` — it must be
`Normal`-ised to a dense array first. `Normal` refuses (leaving the call unevaluated)
above `SA_MAX_ELEMENTS = 2^27` slots or `SA_MAX_RANK = 32`, to avoid a multi-gigabyte
allocation. This is a deliberate, documented limitation.

- Mathilda has **no sparse storage**: a `SparseArray[...]` expression stays
  inert, exactly as written. `Normal` converts it to the dense nested List it
  denotes; that is the conversion Mathematica code relies on, and until
  2026-09-28 `Normal` handed the `SparseArray` back unchanged.
- Positions are Integer lists (an Integer for a vector); without `dims` the
  size is the per-coordinate maximum of the positions. Negative indices count
  from the end once `dims` are known. When several rules name one position the
  **first** wins, as in Mathematica.
- Rule forms: `{p1, p2, ...} -> {v1, v2, ...}`; pattern rules
  (`{i_, i_} -> 1`, `{i_, j_} :> i + j`, conditions allowed), which need `dims`;
  and `Band[start] -> v`, `Band[start, end]`, `Band[start, end, step]`, with `v`
  a scalar or a List of successive values.
- The InputForm of Mathematica's own representation,
  `SparseArray[Automatic, dims, default, {1, {rowptr, colidx}, vals}]`, is
  accepted, so a value pasted out of Mathematica round-trips through `Normal`.
- A specification Mathilda cannot densify (pattern rules without `dims`, a
  position outside `dims`, a result over 2^27 entries) leaves `Normal[...]`
  unevaluated.

**Attributes:** none registered.

## References

**See also:** [Normal](../../data-structures/Normal/)

- Source: [`src/sparsearray.c`](https://github.com/stblake/mathilda/blob/main/src/sparsearray.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_sparsearray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sparsearray.c)

## Notes & additional examples

### Notes

Mathilda has **no sparse storage**: `SparseArray[…]` is an inert symbolic
specification that stays unevaluated at the REPL. To get values out, wrap it in
`Normal`, which builds the dense nested `List` the specification denotes. Because the
object never becomes a compressed buffer, arithmetic and linear algebra are not
accelerated on it — convert with `Normal` first. `Normal` declines (and leaves the
call unevaluated) when the dense array would exceed 2^27 elements or rank 32.
