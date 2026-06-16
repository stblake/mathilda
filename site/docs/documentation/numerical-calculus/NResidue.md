# NResidue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NResidue[expr, {z, z0}]
    numerically finds the residue of expr near z = z0 (the coefficient of (z - z0)^-1 in the Laurent expansion) by integrating around a small circle in the complex plane.
NResidue[{e1, e2, ...}, {z, z0}]
    threads element-wise over the first argument.

Works for essential singularities where the symbolic Residue (which needs a power series) cannot. Cannot distinguish a tiny spurious residual from a true zero -- Chop the result when needed; returns an incorrect value if the contour encloses another singularity or crosses a branch cut.

Options: Radius (contour radius, default 1/100, or Automatic), WorkingPrecision, PrecisionGoal, MaxRecursion (max contour refinements, default 10), Method ('Trapezoidal').
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NResidue[1/x, {x, 0}]
Out[1]= 1.0 + 1.70484e-17*I

In[2]:= NResidue[Sin[1/(10 x)], {x, 0}] // Chop
Out[2]= 0.1

In[3]:= NResidue[1/(1.7 - 2.7 z + z^2), {z, 1.}] // Chop
Out[3]= -1.42857

In[4]:= NResidue[Exp[1/x], {x, 0}, Radius -> 1] // Chop
Out[4]= 1.0

In[5]:= NResidue[{Exp[1/x], Sin[1/x], Cos[1/x]}, {x, 0}, Radius -> 1] // Chop
Out[5]= {1.0, 1.0, 0}

In[6]:= NResidue[1/x + 1/(x + 0.005), {x, 0}, Radius -> 0.001] // Chop
Out[6]= 1.0

In[7]:= NResidue[Exp[1/x], {x, 0}, Radius -> Automatic] // Chop
Out[7]= 1.0

In[8]:= 10! NResidue[Zeta[x]/x^11, {x, 0}, Radius -> 1/2, WorkingPrecision -> 30]
Out[8]= -3628799.999456765884220291526686 + 2.823257449139100034435257009205e-22*I
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)

## Notes & additional examples

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
