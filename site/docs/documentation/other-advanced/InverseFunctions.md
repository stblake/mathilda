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

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Solve[Sin[x] == 1/2, x]
Out[1]= {{x -> ConditionalExpression[2 C[1] Pi + 5/6 Pi, Element[C[1], Integers]]}, {x -> ConditionalExpression[2 C[1] Pi + 1/6 Pi, Element[C[1], Integers]]}}

In[2]:= Solve[Sin[x] == 1/2, x, InverseFunctions -> False]
Out[2]= Solve[Sin[x] == 1/2, x, InverseFunctions -> False]
```

### Notes

`InverseFunctions` is a `Solve` option (default `Automatic`, i.e. enabled) that
controls the inverse-function specialist for invertible heads such as `Sin`,
`Log`, `Exp`, and integer `Power`. With it on, `Sin[x] == 1/2` is inverted to
the full periodic solution set; setting `InverseFunctions -> False` disables the
specialist, so the equation is returned unevaluated.
