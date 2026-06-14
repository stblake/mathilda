# LoadModule

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
LoadModule["relpath"]
    loads the internal Mathilda source module at relpath (relative to
    src/internal), resolving the location independently of the current
    working directory. Each module is loaded at most once. Returns True
    if the module was located and loaded, False otherwise.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`. Returns `True` if the module was located and loaded (or had already

**Attributes:** `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/loadmodule.c`](https://github.com/stblake/mathilda/blob/main/src/loadmodule.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
