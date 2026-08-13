# Residue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Residue[f, {z, z0}]`**

gives the residue of f at the isolated singularity z = z0 -- the coefficient of (z - z0)^-1 in the Laurent expansion of f.

<details>
<summary>Notes</summary>

Computed by power-series expansion, so a residue is found only where f admits a Laurent series at z0. Returns unevaluated at branch points (fractional-power expansions) and when no series can be produced. See NResidue for a numerical alternative that also handles essential singularities.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= Residue[1/z, {z, 0}]
Out[1]= 1

In[2]:= Residue[1/z^2, {z, 0}]
Out[2]= 0

In[3]:= Residue[1/Sin[z]^5, {z, 0}]
Out[3]= 3/8
```

Order-2 pole

```mathematica
In[4]:= Residue[(z + 1)/(z^2 (z - 2)), {z, 0}]
Out[4]= -3/4
```

Complex pole

```mathematica
In[5]:= Residue[1/(z^2 + 1), {z, I}]
Out[5]= -1/2*I
```

Algebraic pole

```mathematica
In[6]:= Residue[x^3/(x^4 - 2), {x, 2^(1/4)}]
Out[6]= 1/4
```

Unknown numerator

```mathematica
In[7]:= Residue[f[z]/z^5, {z, 0}]
Out[7]= 1/24 Derivative[4][f][0]
```

```mathematica
In[8]:= Residue[Zeta[z]/(z - 1)^10, {z, 1}]
Out[8]= -1/362880 StieltjesGamma[9]
```

## Algorithm

residue.c -- Residue[expr, {z, z0}], the symbolic residue.

The residue of f at an isolated singularity z = z0 is the coefficient of (z - z0)^-1 in the Laurent expansion of f. We obtain it directly from the series engine: expand f to order (z - z0)^0 (which always spans the -1 term, however deep the pole), then read the coefficient at exponent -1 out of the resulting SeriesData[z, z0, {coefs}, nmin, nmax, den].

A residue is well defined only for an ordinary Laurent expansion (den == 1). A fractional-power (Puiseux) expansion, den > 1, signals a branch point, where the residue is undefined -- we leave the call unevaluated, matching Mathematica (e.g. Residue[1/Sqrt[z], {z, 0}]).

```text
Algebraic pole locations.  The series engine decides whether z0 is a pole by
```

evaluating the denominator there and testing it against zero; but for a pole whose location is a SUM of radicals (e.g. z0 = -2 + Sqrt[3], a root of 1 + 4 z + z^2), Denominator(z0) is an expression like 1 + 4 (-2 + Sqrt[3]) + (-2 + Sqrt[3])^2 that does not auto-simplify to 0, so

```text
the pole is missed and the residue wrongly comes out 0.  We defeat this by
```

expanding about z0 EXPLICITLY: substitute z -> z0 + w, then Expand the denominator of the result -- polynomial expansion collapses the radical arithmetic (Sqrt[3]^2 -> 3, ...) so the vanishing constant term becomes a

```text
literal 0 and the w-factor of the pole is exposed.  Reading the (z-z0)^-1
```

coefficient is then a plain Series-at-0 of the expanded form.

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [NResidue](../../numerical-calculus/NResidue/), [Together](../../algebra/Together/), [Zeta](../../special-functions/Zeta/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_residue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_residue.c)
