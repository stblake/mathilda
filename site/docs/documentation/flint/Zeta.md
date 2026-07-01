# FLINT`Zeta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FLINT`Zeta[s] gives the numeric value of the Riemann zeta function at the numeric argument s (real or complex), computed to the precision of s (machine precision for exact s) via FLINT's rigorous acb arithmetic (acb_dirichlet_zeta). Unevaluated for symbolic s or at the pole s = 1.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FLINT`Zeta[2]
Out[1]= 1.644934066848226

In[2]:= FLINT`Zeta[4]
Out[2]= 1.082323233711138
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/flint_num_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/flint_num_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
