# Residue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Residue[f, {z, z0}]
    gives the residue of f at the isolated singularity z = z0 -- the coefficient of (z - z0)^-1 in the Laurent expansion of f.

Computed by power-series expansion, so a residue is found only where f admits a Laurent series at z0. Returns unevaluated at branch points (fractional-power expansions) and when no series can be produced. See NResidue for a numerical alternative that also handles essential singularities.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Residue[1/z, {z, 0}]
Out[1]= 1

In[2]:= Residue[1/z^2, {z, 0}]
Out[2]= 0

In[3]:= Residue[1/Sin[z]^5, {z, 0}]
Out[3]= 3/8

In[4]:= Residue[(z + 1)/(z^2 (z - 2)), {z, 0}]      (* order-2 pole *)
Out[4]= -3/4

In[5]:= Residue[1/(z^2 + 1), {z, I}]                (* complex pole *)
Out[5]= -1/2*I

In[6]:= Residue[x^3/(x^4 - 2), {x, 2^(1/4)}]        (* algebraic pole *)
Out[6]= 1/4

In[7]:= Residue[f[z]/z^5, {z, 0}]                   (* unknown numerator *)
Out[7]= 1/24 Derivative[4][f][0]

In[8]:= Residue[Zeta[z]/(z - 1)^10, {z, 1}]
Out[8]= -1/362880 StieltjesGamma[9]
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
