# InverseFunctions

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
InverseFunctions is an option for Solve that enables the
    inverse-function specialist for elementary invertible heads
    (Log, Exp, Sin, Cos, Tan, ArcSin, ArcCos, Sinh, ..., and
    integer Power).  Default: Automatic (enabled).  Setting it
    to False disables the specialist; equations that can only
    be solved through inversion then return unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`InverseFunctions` is an option symbol for `Solve`, not a callable function. It is interned in `sym_names.c` and documented in `solve.c`; it has no builtin handler. `solve.c`'s option parser recognises `InverseFunctions -> spec`: the default `Automatic` (and any value other than `False`) leaves the `enabled` flag set in the solver state (`solveinv.h`), while `InverseFunctions -> False` disables it. When enabled, `Solve` is allowed to "peel" invertible heads (Exp, Log, trig, powers) off an equation by applying their inverse function — the work is done by the inverse-function specialist registered as `` Solve`SolveInverseFunctions `` in `solveinv.c`. The symbol is consumed only as a configuration key.

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/solve.c`](https://github.com/stblake/mathilda/blob/main/src/solve.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
