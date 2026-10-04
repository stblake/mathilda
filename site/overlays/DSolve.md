### Worked examples

```mathematica
In[1]:= DSolve[y'[x] == y[x], y[x], x]  (* a first-order linear ODE *)
```

```mathematica
In[1]:= DSolve[y'[x] == x y[x], y[x], x]  (* separable, with the general constant C[1] *)
```

```mathematica
In[1]:= DSolve[y''[x] + y[x] == 0, y[x], x]  (* constant-coefficient linear: the harmonic oscillator *)
```

```mathematica
In[1]:= DSolve[y''[x] - y[x] == 0, y[x], x]  (* real exponential fundamental set *)
```

```mathematica
In[1]:= DSolve[{y'[x] == y[x], y[0] == 1}, y[x], x]  (* an initial-value problem fits C[1] *)
```

### Notes

`DSolve[eqn, y[x], x]` returns `{{y[x] -> expr}}` with the solution as an
expression in `x`; `DSolve[eqn, y, x]` instead returns the solution as a pure
`Function`. Arbitrary constants are generated as `C[1]`, `C[2]`, ... (rename them
with `GeneratedParameters`). Initial or boundary conditions supplied as equations
at points (e.g. `y[0] == 1`) are fitted, eliminating the constants.

Like `Integrate`, `DSolve` is a cascade polyalgorithm: it tries a sequence of
methods — ordered so that specific, cheap classifiers (factorable, separable,
linear, Bernoulli, exact, constant-coefficient) run before the heavier
second-order and symmetry machinery (Kovacic, Lie point symmetry) and the
Frobenius / power-series fallbacks. A pinned method is available as
`Method -> "<name>"` (e.g. `"Separable"`, `"Kovacic"`), which dispatches directly
with no fallback. Every returned branch is verified by back-substitution before
it is kept, so a decline yields an unevaluated `DSolve[...]` rather than a wrong
answer. Systems `DSolve[{eqns}, {y1, ...}, x]` and partial differential equations
`DSolve[eqn, u, {x, y}]` are handled by their own branches of the dispatcher.
