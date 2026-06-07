# FullForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FullForm[expr]
    prints expr as its raw internal tree (heads written before arguments
    in functional form, no operator or infix sugar).
FullForm is a wrapper recognised by Print/Out; when an input evaluates
to FullForm[expr] the wrapper is consumed by the printer and does not
appear in the output.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
