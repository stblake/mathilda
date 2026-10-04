### Worked examples

```mathematica
In[1]:= Residue[1/z, {z, 0}]  (* a simple pole at the origin *)
```

```mathematica
In[1]:= Residue[1/z^2, {z, 0}]  (* a double pole has no (z)^-1 term, so residue 0 *)
```

```mathematica
In[1]:= Residue[Cot[z], {z, 0}]  (* a transcendental integrand expanded directly *)
```

```mathematica
In[1]:= Residue[1/(z^2 + 1), {z, I}]  (* a simple pole at a complex location *)
```

```mathematica
In[1]:= Residue[Zeta[s], {s, 1}]  (* residue 1 at the simple pole of Zeta *)
```

```mathematica
In[1]:= Residue[1/Sqrt[z], {z, 0}]  (* a branch point -- undefined, left unevaluated *)
```

### Notes

`Residue[f, {z, z0}]` is the coefficient of `(z - z0)^-1` in the Laurent
expansion of `f` at `z0`. It is computed by series expansion: an analytic point
gives `0`, and the expansion order is raised adaptively until the `-1`
coefficient is resolved, so poles of any order are handled.

For a rational integrand with a simple pole the residue is read off as
`P(z0)/Q'(z0)` without inverting a series, which also keeps an algebraic pole
location (a `z0` that is a nested radical) tractable. Transcendental integrands
(`Cot`, `Zeta` near its pole, an unknown `f[z]/z^n`) are expanded directly about
`z0` so the series engine can use its knowledge of their Laurent series there. A
residue is defined only for an ordinary Laurent expansion: a fractional-power
(Puiseux) expansion signals a branch point, where the residue is undefined and
the call is left unevaluated — matching Mathematica. A numerical companion,
`NResidue`, is available for machine-precision work.
