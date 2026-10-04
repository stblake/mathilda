### Worked examples

```mathematica
In[1]:= Resolve[Exists[x, x^2 == 2]]  (* there is a real square root of 2, so the sentence is True *)
```

```mathematica
In[1]:= Resolve[Exists[x, x^2 + 1 == 0]]  (* no real x solves it, so False *)
```

```mathematica
In[1]:= Resolve[Exists[x, a x == 1], Reals]  (* solvable for x exactly when a != 0 *)
```

```mathematica
In[1]:= Reduce[Exists[y, x == y^2], x]  (* the image of squaring is the nonnegative half-line *)
```

### Notes

`Exists[x, expr]` is the quantified statement that *some* value of `x` makes
`expr` true; `Exists[{x1, x2, ...}, expr]` binds several variables and
`Exists[x, cond, expr]` restricts the witness to values satisfying `cond`.

`Exists` is **inert on its own** — it carries `HoldAll`, keeps its bound
variables symbolic, and does not evaluate until `Reduce` or `Resolve` eliminates
it. Evaluate it with `Resolve[Exists[...]]`, or place it inside `Reduce[...]` to
get the condition on the remaining free variables. Quantified problems are solved
over the **Reals**; an explicit non-`Reals` domain declines.

When every variable is bound (a fully quantified *sentence*), the result is
`True` or `False`. When free variables remain — as in `Reduce[Exists[y, x ==
y^2], x]` — the bound variable `y` is projected away by a Cylindrical Algebraic
Decomposition and the answer is a condition on the free ones (here `x >= 0`).
Alternating quantifier prefixes such as `Exists[x, ForAll[y, ...]]` are
eliminated inner-block-first. An undecidable sign or an unsupported construct
leaves the statement unevaluated rather than guessed.
