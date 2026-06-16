### Worked examples

A piecewise definition stays symbolic until its conditions can be decided —
here the absolute-value function, which holds in unevaluated form for symbolic
`x`:

```mathematica
In[1]:= Piecewise[{{-x, x < 0}, {x, x >= 0}}]
Out[1]= Piecewise[{{-x, x < 0}, {x, x >= 0}}, 0]
```

Supplying a concrete value selects the matching branch (and falls through to the
default when no clause applies):

```mathematica
In[1]:= Piecewise[{{1, x > 0}, {-1, x < 0}}, 0] /. x -> 5
Out[1]= 1

In[2]:= Piecewise[{{1, x > 0}}] /. x -> -2
Out[2]= 0
```

Clauses with structurally identical values are automatically merged, their
conditions combined with `Or` — a genuine canonicalisation, not just storage:

```mathematica
In[1]:= Piecewise[{{a, x == 1}, {a, x == 2}, {b, x == 3}}]
Out[1]= Piecewise[{{a, x == 1 || x == 2}, {b, x == 3}}, 0]
```

`Piecewise` integrates with calculus: differentiation threads through every
branch, returning a new piecewise function of the derivatives:

```mathematica
In[1]:= D[Piecewise[{{x^2, x < 0}, {x^3, x >= 0}}], x]
Out[1]= Piecewise[{{2 x, x < 0}, {3 x^2, x >= 0}}, 0]
```

### Notes

`Piecewise[{{v1, c1}, {v2, c2}, ...}, def]` is the symbolic conditional: the `ci`
are tested in order and the value of the first `True` condition is returned, or
`def` (default `0`) if all are `False`. If any earlier condition is not literally
`False` the whole expression is held symbolic. It has attribute `HoldAll`, so only
the `vi` actually returned are evaluated. On construction it canonicalises:
`{vi, False}` clauses are dropped, everything after the first `{vi, True}` is
discarded, and consecutive clauses with identical values are merged via `Or`.
Differentiation and the elementary symbolic machinery thread through the branches.
