### Worked examples

```mathematica
In[1]:= NMaximize[-x^2 + 4 x, x]  (* a downward parabola, maximum at the vertex *)
Out[1]= {4.0, {x -> 2.0}}
```

```mathematica
In[1]:= NMaximize[4 - (x - 1)^2 - (y - 2)^2, {x, y}]  (* a concave bowl, maximum at its centre *)
Out[1]= {4.0, {x -> 1.0, y -> 2.0}}
```

```mathematica
In[1]:= NMaximize[{x y, x + y == 10}, {x, y}]  (* the arithmetic-geometric-mean extremum on x + y = 10 *)
Out[1]= {25.0, {x -> 5.0, y -> 5.0}}
```

### Notes

`NMaximize[f, vars]` searches for a global maximum and returns
`{fmax, {x -> xmax, ...}}`. It is implemented by **minimising `-f` and negating
the objective value**, so it shares `NMinimize`'s methods, options, and
constraint/domain handling in full — see the `NMinimize` page for the method
catalogue (`DifferentialEvolution` by default), the variable and constraint
grammar, integer domains, auto-compilation at `MachinePrecision`, and the
fixed-seed determinism these results rely on. The constrained example attains
`x y = 25` at `x = y = 5`, the classic AM-GM extremum on `x + y = 10`.
`NMaximize` is `Protected` but not `HoldAll`.
