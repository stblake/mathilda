# Dt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Dt[f] gives the total derivative of f.
Dt[f, x] gives the total derivative of f with respect to x.
Dt[f, {x, n}] gives the nth total derivative.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Dt[y^2 + Sin[x]]
Out[1]= Cos[x] Dt[x] + 2 y Dt[y]

In[2]:= Dt[Pi + 3 + x y]
Out[2]= Dt[x] y + x Dt[y]

In[3]:= Dt[y^2, x]
Out[3]= 0

In[4]:= Dt[x^2, {x, 2}]
Out[4]= 2
```

## Implementation notes

**Algorithm.** `Dt` shares the native differentiation core with `D`
(`compute_deriv` in src/calculus/deriv.c). `builtin_dt` handles two modes.
`Dt[f]` (one argument) computes the **total derivative**: it calls
`compute_deriv(f, NULL, NULL)` with a NULL differentiation variable, so unknown
symbols are *not* treated as constants — each contributes a `Dt[sym]`
differential term, and the usual product/quotient/chain rules thread through.
`Dt[f, var, ...]` is defined to be identical to `D[f, var, ...]` (the partial
derivative) and is forwarded to the same per-spec loop used by `builtin_d`
(`parse_var_spec` + `higher_order_partial` / `array_higher_order` /
`compute_deriv_symbolic_order`). Malformed specs emit a `D::dvar`-style message
and return unevaluated.

**Data structures.** `Expr*` tree transformation only; results are returned
un-reduced and folded by the outer evaluator.

**Complexity / limits.** Linear per pass in the tree size. The total-derivative
mode distinguishes itself from `D` solely by the NULL variable that disables the
constant short-circuit; everything else (rules, ownership, fixed-point folding)
is shared with `D`.

- `Protected`, `ReadProtected`.
- Shares the elementary-function derivative table with `D`; the

**Attributes:** `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 2.
- Source: [`src/calculus/deriv.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/deriv.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Dt[x y]
Out[1]= Dt[x] y + x Dt[y]
```

```mathematica
In[1]:= Dt[Sin[x]]
Out[1]= Cos[x] Dt[x]
```

```mathematica
In[1]:= Dt[Log[x]]
Out[1]= Dt[x]/x
```

```mathematica
In[1]:= Dt[a x, x]
Out[1]= a
```

### Notes

`Dt[f]` computes the total differential, treating every symbol as a potential independent variable and emitting `Dt[var]` factors for each one — so `Dt[x y]` gives the full product-rule expansion `Dt[x] y + x Dt[y]`. The two-argument form `Dt[f, x]` is the total derivative with respect to `x`, where other symbols are taken as constants unless they implicitly depend on `x`; `Dt[a x, x]` therefore returns `a`. Elementary functions differentiate through the chain rule with a residual `Dt[x]` factor, as in `Dt[Sin[x]]` and `Dt[Log[x]]`. Constants differentiate to `0` (`Dt[c, x]` gives `0`).
