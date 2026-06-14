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
