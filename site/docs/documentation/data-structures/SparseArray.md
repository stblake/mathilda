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

## Examples (4)

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

## Implementation notes

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

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_sparsearray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sparsearray.c)
