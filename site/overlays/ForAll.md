### Worked examples

```mathematica
In[1]:= Resolve[ForAll[x, x^2 >= 0]]  (* every real square is nonnegative, so True *)
```

```mathematica
In[1]:= Resolve[ForAll[x, a x^2 + 1 > 0], Reals]  (* the parabola stays above zero exactly when a >= 0 *)
```

```mathematica
In[1]:= Resolve[ForAll[x, x^2 + a x + 1 > 0], Reals]  (* positive-definite on the discriminant band -2 < a < 2 *)
```

```mathematica
In[1]:= Reduce[ForAll[y, x <= y^2], x, Reals]  (* below every square iff x <= 0 *)
```

### Notes

`ForAll[x, expr]` is the quantified statement that `expr` holds for *every* value
of `x`; `ForAll[{x1, x2, ...}, expr]` binds several variables and
`ForAll[x, cond, expr]` quantifies over values satisfying `cond` (read as
`cond` implies `expr`).

Like `Exists`, `ForAll` is **inert on its own** (`HoldAll`) and is evaluated only
when `Reduce` or `Resolve` eliminates it, over the **Reals**. A fully quantified
sentence returns `True` or `False`; a statement with free parameters returns the
condition on them — the classic use being the parameter band under which a
family of polynomials keeps one sign, as in the positive-definite quadratic
`x^2 + a x + 1 > 0` holding for all `x` exactly on `-2 < a < 2`.

Internally a universal is decided by the dual emptiness question to the
existential: `ForAll[x, g]` is true exactly when `!g` has no solution. Parametric
cases are projected with a Cylindrical Algebraic Decomposition, and alternating
prefixes are eliminated inner-block-first. An undecidable sign or an unsupported
construct leaves the statement unevaluated.
