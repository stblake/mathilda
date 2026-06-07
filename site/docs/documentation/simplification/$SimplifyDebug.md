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

A diagnostic flag (not a builtin), given an OwnValue defaulting to `False` in
`simp_init`. When set to `True`, `simp_debug_enabled` (read directly off the
OwnValue list to avoid re-evaluation) causes `traced_call_unary` /
`simp_debug_log` to emit one stderr line per transform invocation inside the
Simplify search, in the form `/<TransformName>/: <input> -> <output> [<ms> ms]`,
used to diagnose slow Simplify calls and runaway candidate explosion.

- Default `False`. When set to `True`, `Simplify` prints one line per transform

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/simp/simp_util.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp_util.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
