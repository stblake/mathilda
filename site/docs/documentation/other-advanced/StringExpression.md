# StringExpression

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringExpression[p1, p2, ...] or p1 ~~ p2 ~~ ...`**

Represents a sequence of string patterns to be matched consecutively.

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Flat`, `OneIdentity`, `Protected`.

## References

- Source: [`src/strings/regex/regex_init.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/regex_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringposition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringposition.c)
