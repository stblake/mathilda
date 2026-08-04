# CompilePrint

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
CompilePrint[cf] prints the bytecode of the CompiledFunction cf: its argument and result registers with their types, the scalar/array/tile register banks, and one line per instruction giving both the raw operands and a readable rendering. For an object whose body did not compile it reports the bail reason instead. Returns Null.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/compile/compiled_function.c`](https://github.com/stblake/mathilda/blob/main/src/compile/compiled_function.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
