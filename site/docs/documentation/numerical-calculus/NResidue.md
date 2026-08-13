# NResidue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NResidue[expr, {z, z0}]`**

numerically finds the residue of expr near z = z0 (the coefficient of (z - z0)^-1 in the Laurent expansion) by integrating around a small circle in the complex plane.

**`NResidue[{e1, e2, ...}, {z, z0}]`**

threads element-wise over the first argument.

<details>
<summary>Notes</summary>

Works for essential singularities where the symbolic Residue (which needs a power series) cannot. Cannot distinguish a tiny spurious residual from a true zero -- Chop the result when needed; returns an incorrect value if the contour encloses another singularity or crosses a branch cut. Options: Radius (contour radius, default 1/100, or Automatic), WorkingPrecision, AccuracyGoal (default MachinePrecision), PrecisionGoal, MaxRecursion (max contour refinements, default 10), Method ('Trapezoidal').

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= NResidue[1/x, {x, 0}]
Out[1]= 1.0 + 3.25915e-17*I

In[2]:= NResidue[Sin[1/(10 x)], {x, 0}] // Chop
Out[2]= 0.1

In[3]:= NResidue[1/(1.7 - 2.7 z + z^2), {z, 1.}] // Chop
Out[3]= -1.42857
```

### Options (5)

```mathematica
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

### Applications (4)

```mathematica
In[9]:= Chop[NResidue[Cot[z], {z, 0}]]
Out[9]= 1.0

In[10]:= NResidue[1/(z^2 + 1), {z, I}]
Out[10]= 5.10703e-17 - 0.5*I

In[11]:= NResidue[Gamma[z], {z, -3}]
Out[11]= -0.166667 + 1.41553e-17*I

In[12]:= NResidue[Exp[1/z], {z, 0}, Radius -> 1]
Out[12]= 1.0 - 1.38778e-17*I
```

## Algorithm

nresidue.c — NResidue[expr, {z, z0}, opts]

Numerically finds the residue of `expr` at z = z0 by integrating around a small circle in the complex plane (the periodic-trapezoidal Cauchy integral; see quadrature.h). Unlike the symbolic Residue, this works for essential singularities (Exp[1/x], Sin[1/x]) where the power series does not exist.

```text
  NResidue[expr, {z, z0}]               residue near z = z0
  NResidue[{e1, e2, ...}, {z, z0}]      threads element-wise over arg 1
```

Options (trailing Rule[...] in any order):

```text
  Radius           -> r | Automatic   contour radius (default 1/100)
  WorkingPrecision -> MachinePrecision | digits
  PrecisionGoal    -> digits | Automatic | Infinity
  MaxRecursion     -> n               max N-doublings (default 10)
  Method           -> "Trapezoidal"   (only method in this version)
```

The sampler binds the variable to each complex sample point (Block-style, via temporary OwnValues) and numerically evaluates the integrand, so the user's global symbol table is left untouched. NResidue cannot tell a tiny spurious residual from a true zero — Chop the result when needed.

Memory: receives `res` owned by the evaluator. Returns a fresh Expr* on success or NULL (unevaluated). Never frees `res`. All temporary OwnValues are removed before returning, on every path.

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Residue](../../calculus/Residue/), [Chop](../../elementary-functions/Chop/), [N](../../arithmetic/N/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_nderiv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nderiv.c)
- Tests: [`tests/test_nresidue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nresidue.c)
- Tests: [`tests/test_residue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_residue.c)

## Notes & additional examples

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
