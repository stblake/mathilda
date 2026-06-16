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
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Solve[x^4 + x + 1 == 0, x]
Out[1]= {{x -> Root[1 + #1 + #1^4 &, 1]}, {x -> Root[1 + #1 + #1^4 &, 2]}, {x -> Root[1 + #1 + #1^4 &, 3]}, {x -> Root[1 + #1 + #1^4 &, 4]}}

In[2]:= Quartics
Out[2]= Quartics
```

A biquadratic quartic factors automatically into explicit radicals regardless of the option setting:

```mathematica
In[1]:= Solve[x^4 - 5 x^2 + 6 == 0, x]
Out[1]= {{x -> -Sqrt[2]}, {x -> Sqrt[2]}, {x -> -Sqrt[3]}, {x -> Sqrt[3]}}
```

A pure fourth power yields the four complex fourth roots in closed form:

```mathematica
In[1]:= Solve[x^4 - 2 == 0, x]
Out[1]= {{x -> -2^(1/4)}, {x -> 2^(1/4)}, {x -> -I 2^(1/4)}, {x -> I 2^(1/4)}}
```

A non-biquadratic but radical-solvable quartic gives nested radicals; here the roots are the four values `±Sqrt[2] ± Sqrt[3]` written as `Sqrt[(10 ± 4 Sqrt[6])/2]`:

```mathematica
In[1]:= Solve[x^4 - 10 x^2 + 1 == 0, x]
Out[1]= {{x -> -Sqrt[1/2 (10 - 4 Sqrt[6])]}, {x -> Sqrt[1/2 (10 - 4 Sqrt[6])]}, {x -> -Sqrt[1/2 (10 + 4 Sqrt[6])]}, {x -> Sqrt[1/2 (10 + 4 Sqrt[6])]}}
```

### Notes

`Quartics` is a `Solve` option (the quartic analogue of `Cubics`), not a
function; evaluating the bare symbol just returns itself. With the default
`Quartics -> False`, an irreducible quartic is returned as held `Root[]`
objects as above; `Quartics -> True` requests explicit radical formulas where
they apply (biquadratic and other special quartics still reduce automatically).
