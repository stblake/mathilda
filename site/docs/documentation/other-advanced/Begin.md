# Begin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Begin["ctx`"] sets the current context ($Context) to "ctx`", saving
the previous value for End[] to restore. If the argument starts with a
backtick it is interpreted relative to the current context:
Begin["`Private`"] inside "MyPkg`" yields "MyPkg`Private`".
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
