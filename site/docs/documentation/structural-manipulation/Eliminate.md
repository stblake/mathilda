# Eliminate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Eliminate[eqns, vars]
    eliminates vars between a list/conjunction of simultaneous
    equations lhs == rhs, returning a balanced Equal[] or an
    And[] of Equal[]s in the remaining variables (True if the
    elimination ideal is empty, False if the system is
    inconsistent). Works on polynomial equations over Q via a
    lexicographic Gröbner basis with elimination block. A
    principal-branch inverse-function pre-pass peels single Sin/
    Cos/Tan/Sinh/Cosh/Tanh/Exp/Log wrappers and emits
    Eliminate::ifun; non-polynomial systems otherwise return
    unevaluated with Eliminate::nlin.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- Drives the lex-order Buchberger engine (`GroebnerBasis`) with an

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
