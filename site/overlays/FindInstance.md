### Worked examples

```mathematica
In[1]:= FindInstance[x + y == 10 && x - y == 2, {x, y}]  (* one witness to a linear system *)
```

```mathematica
In[1]:= FindInstance[x^2 == 2, x]  (* a single root of a quadratic over the complexes *)
```

```mathematica
In[1]:= FindInstance[x^2 + y^2 == 1 && x > 0 && y > 0, {x, y}, Reals]  (* a point on the unit circle in the open first quadrant *)
```

```mathematica
In[1]:= FindInstance[2 x + 3 y == 1, {x, y}, Integers]  (* an integer lattice point on the line *)
```

```mathematica
In[1]:= FindInstance[x^2 + 1 == 0 && x > 0, x, Reals]  (* provably empty: returns {} *)
```

```mathematica
In[1]:= FindInstance[a && Xor[a, b], {a, b}, Booleans]  (* a satisfying Boolean assignment *)
```

### Notes

`FindInstance[expr, vars]` returns *one* instance of `vars` that satisfies the
statement `expr`, in `Solve`'s rule-list form `{{x -> v, ...}}`; `{}` means the
solution set is provably empty. `FindInstance[expr, vars, dom]` names the domain
(`Complexes`, `Reals`, `Integers`, `Rationals`, or `Booleans`), and a trailing
integer `n` asks for up to `n` instances. The default domain is `Complexes`, or
`Reals` when `expr` carries an ordering — the same rule `Reduce` uses.

**Every instance returned is verified** against the original statement, so a
reported point is always a true solution. Because of that, `FindInstance` can
succeed where `Reduce` gives no complete reduction: it instantiates parametric
Diophantine families, searches a bounded integer box over `Integers`, finds
branch-cut and open-region witnesses by structured sampling, and falls back to a
numerical feasibility search for transcendental or inexact real systems.
Variables may be plain symbols or indexed forms `c[i]`.

A returned `{}` is a genuine emptiness proof — from a decidable `Reduce`, an
exhausted finite integer box, or a Gröbner certificate — whereas a system where
no witness is found and emptiness is not proved stays unevaluated. Use
`Modulus -> p` to search over `Z/pZ`.
