# Quartics

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
Quartics is an option for Solve that controls whether quartic
    equations are solved via explicit radical formulas
    (Quartics -> True) or returned as held Root[] objects
    (default Quartics -> False).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Quartics` is not a builtin function — it is an option symbol for `Solve`, registered (with its docstring) in `solve_init`. It is recognised by `is_known_option_name` and consumed in `apply_option`, which sets `opts.poly.quartics_radical = (rhs === True)`. With the default `Quartics -> False`, quartic equations are returned as held `Root[]` objects; `Quartics -> True` requests explicit radical (Ferrari) formulas from the polynomial specialist (`src/poly/solvepoly.c`). The same option is forwarded by `Eigenvalues`/`Eigensystem` (`src/linalg/`) so the characteristic-polynomial roots can likewise be returned as radicals or held `Root[]`s.

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/solve.c`](https://github.com/stblake/mathilda/blob/main/src/solve.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Solve[x^4 + x + 1 == 0, x]
Out[1]= {{x -> Root[1 + #1 + #1^4 &, 1]}, {x -> Root[1 + #1 + #1^4 &, 2]}, {x -> Root[1 + #1 + #1^4 &, 3]}, {x -> Root[1 + #1 + #1^4 &, 4]}}

In[2]:= Quartics
Out[2]= Quartics
```

### Notes

`Quartics` is a `Solve` option (the quartic analogue of `Cubics`), not a
function; evaluating the bare symbol just returns itself. With the default
`Quartics -> False`, an irreducible quartic is returned as held `Root[]`
objects as above; `Quartics -> True` requests explicit radical formulas where
they apply (biquadratic and other special quartics still reduce automatically).
