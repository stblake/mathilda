# FactorTermsList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FactorTermsList[poly]
    gives a list in which the first element is the overall numerical factor in poly, and the second element is the polynomial with the overall factor removed.
FactorTermsList[poly, {x1, x2, ...}]
    gives a list of factors of poly. The first element in the list is the overall numerical factor. The second element is a factor that does not depend on any of the xi. Subsequent elements are factors which depend on progressively more of the xi.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- The product of the returned list always reproduces the input (after

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
