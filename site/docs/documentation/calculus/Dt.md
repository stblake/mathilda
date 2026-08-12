# Dt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Dt[f] gives the total derivative of f.`**

**`Dt[f, x] gives the total derivative of f with respect to x.`**

**`Dt[f, {x, n}] gives the nth total derivative.`**

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Dt[y^2 + Sin[x]]
Out[1]= Cos[x] Dt[x] + 2 y Dt[y]

In[2]:= Dt[Pi + 3 + x y]
Out[2]= Dt[x] y + x Dt[y]

In[3]:= Dt[x^n, x]
Out[3]= x^(-1 + n) (n + x Log[x] Dt[n, x])

In[4]:= Dt[a x + b, x]
Out[4]= a + x Dt[a, x] + Dt[b, x]

In[5]:= Dt[x^2, {x, 2}]
Out[5]= 2
```

### Applications (7)

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

```mathematica
In[1]:= Dt[x^n]
Out[1]= x^(-1 + n) (n Dt[x] + Dt[n] x Log[x])
```

```mathematica
In[1]:= Dt[f[g[x]]]
Out[1]= Dt[x] Derivative[1][g][x] Derivative[1][f][g[x]]
```

```mathematica
In[1]:= Dt[x^2 y^3]
Out[1]= 2 x Dt[x] y^3 + 3 x^2 y^2 Dt[y]
```

## Options & behaviour

### Examples

## Algorithm

deriv.c -- Native C implementation of Mathematica-style differentiation.

This module replaces the fragile rule-based bootstrap in src/internal/deriv.m with a direct, dispatch-driven implementation.

Overview -------- The two key entry points are the builtins D (partial derivative) and Dt (total derivative). Both ultimately funnel through a single recursive core, ``compute_deriv``, parameterised by an optional differentiation variable. When the variable is non-NULL we compute a partial derivative treating everything else as constant (using a fast FreeQ-style walk to short-circuit constant sub-trees). When the variable is NULL we compute a total derivative -- unknown symbols then participate as ``Dt[sym]`` terms.

Why this is faster than the rule-based implementation ----------------------------------------------------- The old deriv.m relied on ~60 DownValues. Each call to D[f, x] would:

```text
  * scan the DownValues list for D linearly,
  * attempt pattern matching against every rule head (Plus, Times,
    Power, every elementary function, ...),
  * run ``/;`` side-conditions such as FreeQ,
  * perform attempt-evaluate/backtrack cycles in the matcher,
  * recursively re-evaluate the result through the full rule engine.
```

In contrast, this module performs a single head-symbol strcmp dispatch per call, constructs the derivative expression directly, and lets the outer evaluator simplify arithmetic. Crucially, the constant-detection step uses a tailored structural traversal (expr_free_of) that avoids calling out to the generic FreeQ builtin.

Returned expressions -------------------- Every builder below produces plain un-reduced expression trees (e.g. Plus[0, x] or Times[1, x]). The outer Mathilda evaluator runs a full fixed-point reduction on the value we return, so Plus[0, ...], Times[1, ...], and all subsequent chain-rule simplifications fold automatically. This keeps the code readable and avoids duplicating the arithmetic simplifier.

Memory ownership ---------------- Every helper that returns an ``Expr*`` returns a freshly allocated tree owned by the caller. Input expressions are never mutated; sub-expressions that need to be reused are always deep-copied.

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
  only dispatch difference is the base-case handling of symbols
  (free symbols become `Dt[s, x]` factors instead of `0`).

**Attributes:** `Protected`, `ReadProtected`.

## See also

[Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/), [I](../../mathematical-constants/I/), [EulerGamma](../../mathematical-constants/EulerGamma/), [Catalan](../../mathematical-constants/Catalan/), [GoldenRatio](../../mathematical-constants/GoldenRatio/), [Degree](../../mathematical-constants/Degree/), [D](../../calculus/D/)

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 2.
- Source: [`src/calculus/deriv.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/deriv.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_eliminate.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eliminate.c)
- Tests: [`tests/test_powerexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_powerexpand.c)

## Notes & additional examples

### Notes

`Dt[f]` computes the total differential, treating every symbol as a potential independent variable and emitting `Dt[var]` factors for each one — so `Dt[x y]` gives the full product-rule expansion `Dt[x] y + x Dt[y]`. The two-argument form `Dt[f, x]` is the total derivative with respect to `x`, where other symbols are taken as constants unless they implicitly depend on `x`; `Dt[a x, x]` therefore returns `a`. Elementary functions differentiate through the chain rule with a residual `Dt[x]` factor, as in `Dt[Sin[x]]` and `Dt[Log[x]]`. Constants differentiate to `0` (`Dt[c, x]` gives `0`).
