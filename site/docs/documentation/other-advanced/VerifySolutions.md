# VerifySolutions

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
VerifySolutions is an option for Solve that decides whether to
    verify each returned solution by back-substitution.
    Default: Automatic.  Reserved.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`VerifySolutions` is an option *symbol* for `Solve`, not a callable function. In `src/solve.c` it is listed in `is_known_option_name`, so `VerifySolutions -> _` is peeled off the argument list as a valid option (`is_option_arg`) rather than treated as a variable. However `apply_option` currently does nothing with it — the value is parsed and accepted but not yet wired into the polynomial solver (the docstring notes `Default: Automatic. Reserved.`). It exists so user code can pass the option without error.

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/solve.c`](https://github.com/stblake/mathilda/blob/main/src/solve.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
