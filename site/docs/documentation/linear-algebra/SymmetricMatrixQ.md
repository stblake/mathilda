# SymmetricMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SymmetricMatrixQ[m]
    gives True if m is explicitly symmetric (m == Transpose[m]),
    and False otherwise.

Options:
    SameTest  -> Automatic   function used to test equality of entries.
    Tolerance -> Automatic   numeric tolerance for approximate matrices.

With SameTest -> f, entries m[i,j] and m[j,i] are taken to be equal
when f[m[i,j], m[j,i]] gives True.  With Tolerance -> t, entries are
accepted when Abs[m[i,j] - m[j,i]] <= t.  SymmetricMatrixQ uses the
definition m^T == m for both real- and complex-valued matrices, so a
complex symmetric matrix need not be Hermitian.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- Default test is structural via `expr_eq`; the diagonal is exempt

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
