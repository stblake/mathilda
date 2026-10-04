### Worked examples

```mathematica
In[1]:= NMinimize[x^2 - 4 x + 7, x]  (* a convex parabola, minimum at the vertex *)
Out[1]= {3.0, {x -> 2.0}}
```

```mathematica
In[1]:= NMinimize[{x + y, x^2 + y^2 == 1}, {x, y}]  (* an equality constraint: smallest x + y on the unit circle *)
Out[1]= {-1.41421, {x -> -0.707107, y -> -0.707107}}
```

```mathematica
In[1]:= NMinimize[{2 x + 3 y, x >= 1, y >= 1, x + y >= 4}, {x, y}]  (* a small linear program with inequality constraints *)
Out[1]= {9.0, {x -> 3.0, y -> 1.0}}
```

```mathematica
In[1]:= NMinimize[{(x - 3)^2 + 1, Element[x, Integers]}, x]  (* restricting x to the integers *)
Out[1]= {1.0, {x -> 3}}
```

### Notes

`NMinimize[f, vars]` searches for a global minimum and returns
`{fmin, {x -> xmin, ...}}`. Variables may be bare symbols, `{x, lo, hi}` search
intervals, or indexed variables, and the constraint set may mix equalities,
inequalities, chained inequalities, `And`, and `Or`. A scalar integer variable
is declared with `Element[x, Integers]`.

The default `Method -> Automatic` is `"DifferentialEvolution"` (DE/rand/1/bin
with Deb's feasibility rules); `"NelderMead"`, `"RandomSearch"`,
`"SimulatedAnnealing"`, `"SHGO"`, `"DualAnnealing"`, `"DIRECT"`, and
`"BasinHopping"` are also available. The global-search point is polished by an
exact local optimizer, so smooth convex and polynomial problems return clean
values. The search is **deterministic for a fixed `RandomSeed`** (the default
seed is fixed), so these results are reproducible. At `MachinePrecision` the
objective and constraints are auto-compiled to bytecode for the trial-point
loop; `WorkingPrecision -> d` refines a continuous box/unconstrained problem to
`d` digits with MPFR. An empty feasible set returns
`{Infinity, {x -> Indeterminate, ...}}`. `NMinimize` is `Protected` but not
`HoldAll`; its variables should be unbound symbols, and their values are set and
restored `Block`-style during the search.
