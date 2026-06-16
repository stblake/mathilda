### Worked examples

```mathematica
In[1]:= Chop[NResidue[Cot[z], {z, 0}]]
Out[1]= 1.0
```

```mathematica
In[1]:= NResidue[1/(z^2 + 1), {z, I}]
Out[1]= 5.10703e-17 - 0.5*I
```

```mathematica
In[1]:= NResidue[Gamma[z], {z, -3}]
Out[1]= -0.166667 + 1.41553e-17*I
```

```mathematica
In[1]:= NResidue[Exp[1/z], {z, 0}, Radius -> 1]
Out[1]= 1.0 - 1.38778e-17*I
```

### Notes

`NResidue[expr, {z, z0}]` finds the residue (the coefficient of `(z - z0)^-1` in
the Laurent expansion) by integrating around a small circle. `Cot[z]` has a
simple pole at `0` with residue `1`. The rational function gives the expected
`-I/2` at the pole `z = I`. The Gamma function has a pole at every non-positive
integer with residue `(-1)^k / k!`, so at `z = -3` the value is
`-1/6 = -0.166667...`. The most striking case is the essential singularity of
`Exp[1/z]` at the origin: its residue is `1`, which the symbolic `Residue`
(needing a power series) cannot obtain. Use `Radius` to control the contour and
`Chop` to clear spurious imaginary residuals.
