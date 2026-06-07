# $SimplifyDebug

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
$SimplifyDebug
    When set to True, Simplify prints one stderr line per
    transform invocation: /Name/: <input> -> <output> [<ms> ms].
    Defaults to False. Useful for diagnosing slow Simplify calls.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Default `False`. When set to `True`, `Simplify` prints one line per transform

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/simp/simp.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
