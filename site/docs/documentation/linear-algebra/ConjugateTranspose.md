# ConjugateTranspose

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ConjugateTranspose[m]
    Gives the conjugate transpose of m, equivalent to Conjugate[Transpose[m]].
ConjugateTranspose[m, spec]
    Gives Conjugate[Transpose[m, spec]], permuting the levels of m according to
    the spec list and then conjugating every entry.
    On a 1-D vector, ConjugateTranspose[vec] conjugates the entries without
    changing the shape of vec.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ConjugateTranspose[{{1, 2 I, 3}, {3 + 4 I, 5, I}}]
Out[1]= {{1, 3 - 4*I}, {-2*I, 5}, {3, -I}}

In[2]:= ConjugateTranspose[{{a + b I, c + d I}}]
Out[2]= {{Conjugate[a + I b]}, {Conjugate[c + I d]}}

In[3]:= ConjugateTranspose[{1, 2 I, 3 + 4 I}]
Out[3]= {1, -2*I, 3 - 4*I}
```

## Implementation notes

`builtin_conjugate_transpose` (in `src/list.c`) is a thin composition over existing primitives. It first checks the argument is a rectangular nested `List` via `get_array_dimensions`; a symbolic (non-list) matrix is left unevaluated so `ConjugateTranspose[A]` survives. For a 1-D vector it just maps `Conjugate` elementwise. Otherwise it builds and evaluates `Transpose[m]` (or `Transpose[m, spec]`), then conjugates the transposed result. All heavy lifting is delegated to `Transpose` and `Conjugate` through `eval_and_free`.

- `Protected`.
- On a 1-D vector, `ConjugateTranspose[vec]` conjugates the entries but

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
