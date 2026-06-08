# GeneratedParameters

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
GeneratedParameters is an option for Solve specifying the
    head used for fresh integer-parameter symbols introduced by
    the inverse-function specialist.  Default: C, giving
    C[1], C[2], ...  Only the bare-symbol form is honoured;
    the Function form is reserved.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`GeneratedParameters` is an option symbol for `Solve`, not a callable function. It is interned in `sym_names.c` and given a docstring in `solve.c`; it has no builtin handler. `solve.c`'s option parser recognises `GeneratedParameters -> head` and stores the (interned) head string in the solver state's `param_head` field (declared in `solveinv.h`, default `"C"`). When `Solve` returns a parameterised family (e.g. integer-multiple solutions of a trig equation), the inverse-function specialist (`solveinv.c`) emits the free parameters as `head[k]` — so `GeneratedParameters -> C` produces `C[1]`, `C[2]`, … This is purely a configuration key consulted during solving.

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/solve.c`](https://github.com/stblake/mathilda/blob/main/src/solve.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Solve[Exp[x] == 2, x]
Out[1]= {{x -> ConditionalExpression[Log[2] + (2*I) C[1] Pi, Element[C[1], Integers]]}}

In[2]:= Solve[Exp[x] == 2, x, GeneratedParameters -> K]
Out[2]= {{x -> ConditionalExpression[Log[2] + (2*I) K[1] Pi, Element[K[1], Integers]]}}
```

### Notes

`GeneratedParameters` is a `Solve` option that names the head used for the
fresh integer-parameter symbols the inverse-function solver introduces for
equations with infinitely many solutions. The default `C` produces
`C[1], C[2], ...`; here `GeneratedParameters -> K` switches them to `K[1]`.
