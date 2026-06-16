# D

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
D[f, x] gives the partial derivative of f with respect to x.
D[f, {x, n}] gives the nth partial derivative.
D[f, x, y, ...] gives the mixed derivative.
D[f, x, NonConstants -> {y, ...}] treats the listed symbols as implicit functions of x.
Distributes over Equal: D[a == b, x] gives D[a, x] == D[b, x].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= D[x^3, x]
Out[1]= 3 x^2

In[2]:= D[Sin[x^2], x]
Out[2]= 2 x Cos[x^2]

In[3]:= D[Sin[a x], {x, 3}]
Out[3]= -a^3 Cos[a x]

In[4]:= D[f[g[x]], x]
Out[4]= Derivative[1][g][x] Derivative[1][f][g[x]]

In[5]:= D[f[x, y], x]
Out[5]= Derivative[1, 0][f][x, y]

In[6]:= D[Derivative[2][f][x], x]
Out[6]= Derivative[3][f][x]

In[7]:= D[{x, x^2, Sin[x]}, x]
Out[7]= {1, 2 x, Cos[x]}

In[8]:= D[Log[b, x], x]
Out[8]= 1/(Log[b] x)
```

## Implementation notes

**Algorithm.** `D[f, x]` is computed by a native, dispatch-driven C
differentiator (`builtin_d` → `compute_deriv` in src/calculus/deriv.c), which
replaced the old rule-based `src/internal/deriv.m` (now a no-op stub kept only
so users can drop in custom `Dt` identities). `builtin_d` first splits trailing
arguments into `NonConstants -> ...` options and variable specs, then applies
each spec sequentially so that mixed partials `D[f, x, y]` and higher orders
`D[f, {x, n}]` and array forms `D[f, {{x1,...,xN}}]` all reduce to repeated
single-variable differentiation (`higher_order_partial`, `array_higher_order`,
and `compute_deriv_symbolic_order` for symbolic order `n`).

The core `compute_deriv(f, x, nonconsts)` performs a single head-symbol dispatch
per node: linearity for `Plus`; the general n-factor **product rule** for
`Times`; `Power` handled with the combined power/exponential/general
`u^v` rule; the **quotient** case falls out of `Power[..,-1]`; and elementary
unary functions (Sin, Cos, Exp, Log, ArcTan, …) get their derivative from a
table lookup (`elementary_fprime`) composed with the **chain rule**. Unknown
single- and multi-argument functions get the generic chain rule
(`chain_rule_unknown`), emitting `Derivative[...][f]` factors. Constant
subtrees are short-circuited by a tailored structural walk (`expr_free_of`)
rather than calling the generic `FreeQ` builtin — this is the main speedup over
the rule-based version. Each builder returns plain un-reduced trees (e.g.
`Plus[0, x]`, `Times[1, x]`); the outer Mathilda fixed-point evaluator folds the
arithmetic.

**Data structures.** Pure `Expr*` tree transformation; no auxiliary numeric
representation. Symbolic-order derivatives produce closed forms where possible
and otherwise fall back to an unevaluated `D[...]`.

**Complexity / limits.** Linear in the size of the input tree per
differentiation pass (with constant-subtree pruning); `D[f, x, n]` costs n
passes. `NonConstants` is honoured for ordinary specs but not threaded through
the closed-form symbolic-order path.

- `Protected`, `ReadProtected`.
- Recognises the elementary heads `Plus`, `Times`, `Power`, `Sqrt`,

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
In[1]:= D[x^n, x]
Out[1]= n x^(-1 + n)
```

```mathematica
In[1]:= D[Exp[a x], x]
Out[1]= a E^(a x)
```

```mathematica
In[1]:= D[x^2 y, {x, 2}]
Out[1]= 2 y
```

```mathematica
In[1]:= D[x^x, x]
Out[1]= x^(-1 + x) (x + x Log[x])
```

```mathematica
In[1]:= D[Sin[x]^Cos[x], x]
Out[1]= Sin[x]^(-1 + Cos[x]) (Cos[x]^2 - Sin[x]^2 Log[Sin[x]])
```

```mathematica
In[1]:= D[Log[Gamma[x]], x]
Out[1]= PolyGamma[0, x]
```

```mathematica
In[1]:= D[f[g[x]], x]
Out[1]= Derivative[1][g][x] Derivative[1][f][g[x]]
```

### Notes

`D` is driven by the pattern-based derivative rules in `src/internal/deriv.m`, implementing the chain, product, and quotient rules together with the elementary-function table. Logarithmic differentiation is handled automatically, so `D[x^x, x]` and the more general `D[Sin[x]^Cos[x], x]` (variable base *and* variable exponent) come out correctly, and `D[Log[Gamma[x]], x]` is recognised as the digamma function `PolyGamma[0, x]`. Differentiating an unknown function head produces a `Derivative[n][f]` operator rather than evaluating further, so the chain rule on `f[g[x]]` returns a product of such operators. The `{x, n}` form takes the `n`th derivative and treats other symbols as constants by default; use `NonConstants -> {...}` to mark implicit dependencies.
