# FLINT`HurwitzZeta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FLINT`HurwitzZeta[s, a] gives the numeric value of the Hurwitz zeta function via FLINT (acb_dirichlet_hurwitz), to the precision of the arguments. Unevaluated for symbolic arguments or at a pole.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FLINT`HurwitzZeta[2, 1/2]
Out[1]= 4.934802200544679
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/flint_num_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/flint_num_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
