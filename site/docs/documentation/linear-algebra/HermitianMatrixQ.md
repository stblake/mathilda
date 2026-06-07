# HermitianMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HermitianMatrixQ[m]
    gives True if m is explicitly Hermitian (m == ConjugateTranspose[m]),
    and False otherwise.

Options:
    SameTest  -> Automatic   function used to test equality of entries.
    Tolerance -> Automatic   numeric tolerance for approximate matrices.

With SameTest -> f, entries m[i,j] and Conjugate[m[j,i]] are taken to be
equal when f[m[i,j], Conjugate[m[j,i]]] gives True.  With Tolerance -> t,
entries are accepted when Abs[m[i,j] - Conjugate[m[j,i]]] <= t.
Diagonal entries must satisfy the same test (i.e. be purely real for
numeric matrices).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- Default test is structural: it accepts (Conjugate[a], a) / (a, Conjugate[a])

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
