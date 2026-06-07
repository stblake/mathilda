# $PreRead

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$PreRead
    is a global variable whose value, if set, is applied to the
    text of every input expression before it is fed to Mathilda.

The value of $PreRead must be a function that takes one string
argument and returns a string. If unset, the input is parsed
without modification. If the hook returns a non-string value
the original input is used and a $PreRead::strret diagnostic
is emitted.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/repl_hooks.c`](https://github.com/stblake/mathilda/blob/main/src/repl_hooks.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
