### Worked examples

```mathematica
In[1]:= Resolve[Exists[x, a x == 1], Reals]  (* an inverse for a exists exactly when a != 0 *)
```

```mathematica
In[1]:= Resolve[ForAll[x, x^2 + a x + 1 > 0], Reals]  (* discriminant band for a positive-definite quadratic *)
```

```mathematica
In[1]:= Resolve[Exists[x, x^2 == 2 && x > 0]]  (* fully quantified sentence, decides to True *)
```

```mathematica
In[1]:= Resolve[ForAll[x, Exists[y, y > x]]]  (* an alternating prefix, eliminated inner-block-first *)
```

### Notes

`Resolve[expr]` eliminates the quantifiers (`Exists`, `ForAll`) from `expr`,
returning an equivalent quantifier-free statement: `True` or `False` for a fully
quantified sentence, or a condition on the remaining free variables otherwise. It
is the dedicated quantifier-elimination front-end — its argument must be a
quantified head, and it shares the engine `Reduce` uses when its own input happens
to be quantified. `Resolve[expr, dom]` names the domain, but **Reals** is the
default and only supported one.

Alternating quantifier prefixes such as `ForAll[x, Exists[y, y > x]]` are handled
by eliminating the innermost block to a quantifier-free formula and then
re-eliminating the enclosing block over it, composing to arbitrary depth. The
parametric case is driven by a Cylindrical Algebraic Decomposition that projects
away the bound variables and emits a sign formula over the free ones.

`Resolve` follows the engine's soundness rule: an undecidable sign, a non-`Reals`
domain, or a construct outside the supported fragment leaves the call unevaluated
rather than returning a wrong formula.
