# ConjugateTranspose

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ConjugateTranspose[m]`**

Gives the conjugate transpose of m, equivalent to Conjugate\[Transpose\[m\]\].

**`ConjugateTranspose[m, spec]`**

Gives Conjugate\[Transpose\[m, spec\]\], permuting the levels of m according to the spec list and then conjugating every entry. On a 1-D vector, ConjugateTranspose\[vec\] conjugates the entries without changing the shape of vec.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ConjugateTranspose[{{1, 2 I, 3}, {3 + 4 I, 5, I}}]
Out[1]= {{1, 3 - 4*I}, {-2*I, 5}, {3, -I}}

In[2]:= ConjugateTranspose[{{a + b I, c + d I}}]
Out[2]= {{Conjugate[a + I b]}, {Conjugate[c + I d]}}

In[3]:= ConjugateTranspose[{1, 2 I, 3 + 4 I}]
Out[3]= {1, -2*I, 3 - 4*I}
```

### Applications (3)

```mathematica
In[4]:= ConjugateTranspose[{{1 + I, 2}, {3, 4 - I}}]
Out[4]= {{1 - I, 3}, {2, 4 + I}}

In[5]:= ConjugateTranspose[{{a, b}, {c, d}}]
Out[5]= {{Conjugate[a], Conjugate[c]}, {Conjugate[b], Conjugate[d]}}

In[6]:= m = {{1, I}, {-I, 2}}; ConjugateTranspose[m] == m
Out[6]= True
```

## Implementation notes

`builtin_conjugate_transpose` (in `src/list.c`) is a thin composition over existing primitives. It first checks the argument is a rectangular nested `List` via `get_array_dimensions`; a symbolic (non-list) matrix is left unevaluated so `ConjugateTranspose[A]` survives. For a 1-D vector it just maps `Conjugate` elementwise. Otherwise it builds and evaluates `Transpose[m]` (or `Transpose[m, spec]`), then conjugates the transposed result. All heavy lifting is delegated to `Transpose` and `Conjugate` through `eval_and_free`.

- `Protected`.
- On a 1-D vector, `ConjugateTranspose[vec]` conjugates the entries but
  does not change the shape of `vec` (matches the Mathematica convention).
- For symbolic entries, `Conjugate[x]` is left wrapped around `x`.
- Works on higher-rank tensors with the same `spec` semantics as
  `Transpose`.
- Reads and returns a [packed list](../packed-arrays/index.md) at every rank, including
  the rank-1 form. A **real** buffer skips the conjugation entirely — it is the
  identity there, and running it would be a second full pass over the data.

**Attributes:** `Protected`.

## References

**See also:** [Transpose](../../structural-manipulation/Transpose/)

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_compile_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_linalg.c)
- Tests: [`tests/test_conjugate_transpose.c`](https://github.com/stblake/mathilda/blob/main/tests/test_conjugate_transpose.c)
- Tests: [`tests/test_hermitian_matrix_q.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hermitian_matrix_q.c)
- Tests: [`tests/test_negative_definite_matrix_q.c`](https://github.com/stblake/mathilda/blob/main/tests/test_negative_definite_matrix_q.c)

## Notes & additional examples

### Notes

`ConjugateTranspose[m]` is the Hermitian adjoint `Conjugate[Transpose[m]]`. The last example confirms a matrix is Hermitian (equal to its own conjugate transpose). On a vector it conjugates the entries in place.
