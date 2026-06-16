# Eliminate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Eliminate[eqns, vars]
    eliminates vars between a list/conjunction of simultaneous
    equations lhs == rhs, returning a balanced Equal[] or an
    And[] of Equal[]s in the remaining variables (True if the
    elimination ideal is empty, False if the system is
    inconsistent). Works on polynomial equations over Q via a
    lexicographic Gröbner basis with elimination block. A
    principal-branch inverse-function pre-pass peels single Sin/
    Cos/Tan/Sinh/Cosh/Tanh/Exp/Log wrappers and emits
    Eliminate::ifun; non-polynomial systems otherwise return
    unevaluated with Eliminate::nlin.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_eliminate` (in `src/poly/eliminate.c`) removes a set of variables from a system of equations by **Gröbner elimination**. It accepts `Eliminate[eqns, vars]` where `eqns` is a `List`/`And` of `lhs == rhs` equations. A pre-pass tries to handle simple invertible transcendental equations of the shape `f[poly] == const` (e.g. `Exp`, `Log`, trig) by a one-layer principal-branch rewrite, emitting an `Eliminate::ifun` diagnostic to warn that branches may be lost.

The algebraic path moves each equation to `lhs − rhs` form, collects all variables, and orders them so the variables to be eliminated form the leading block (the same elimination-block layout used by `GroebnerBasis`'s 3-arg form). It converts the polynomials to `GBPoly` via `gb_from_expr` and runs `gb_buchberger` under the lex/elimination order. By the elimination theorem, the basis polynomials free of the eliminated variables generate the elimination ideal; those survivors are re-presented as balanced `Equal[posPart, negPart]` equations and joined with `And` when there are several.

**Data structures.** Reuses the Gröbner subsystem's `GBPoly` (GMP `mpq_t` coefficients + row-major exponent matrix) and the Buchberger/Gebauer–Möller core in `groebner.c`. The memory contract is the standard builtin one — every early-return path frees temporaries and never frees `res`.

**Complexity / limits.** Inherits Buchberger's worst-case doubly-exponential cost (lex/elimination orders are the expensive case). `gb_from_expr` cannot atomise Power-headed main variables, so genuinely transcendental systems outside the simple-invertible pre-pass fall back to leaving the input unevaluated.

- `Protected`.
- Drives the lex-order Buchberger engine (`GroebnerBasis`) with an

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- T. Becker, V. Weispfenning, *Gröbner Bases* (Springer, 1993).
- D. Cox, J. Little, D. O'Shea, *Ideals, Varieties, and Algorithms* (Springer).
- Source: [`src/poly/eliminate.c`](https://github.com/stblake/mathilda/blob/main/src/poly/eliminate.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Eliminate[{x + y == 2, x - y == 0}, y]
Out[1]= x == 1

In[2]:= Eliminate[{a == b + c, d == a - c}, c]
Out[2]= d == b
```

```mathematica
In[1]:= Eliminate[{a == x + y, b == x y}, {x, y}]
Out[1]= True

In[2]:= Eliminate[{p == x + 1/x, q == x^2 + 1/x^2}, x]
Out[2]= 2 + q == p^2
```

```mathematica
In[1]:= Eliminate[{x == a Cos[t], y == a Sin[t]}, t]
Out[1]= x^2 + y^2 == a^2
```

```mathematica
In[1]:= Eliminate[{u == Exp[x], v == Exp[2 x]}, x]
Out[1]= v == u^2
```

### Notes

`Eliminate[eqns, vars]` removes `vars` from a system of polynomial equations
over the rationals (via a lexicographic Gröbner basis with an elimination
block), returning the relations among the remaining variables. It yields
`True` if the elimination ideal is empty and `False` if the system is
inconsistent; non-polynomial systems return unevaluated with `Eliminate::nlin`.
