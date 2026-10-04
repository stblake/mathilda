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

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (9)

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

Branch point

```mathematica
In[9]:= Residue[1/Sqrt[z], {z, 0}]
Out[9]= Residue[1/Sqrt[z], {z, 0}]
```

### Applications (6)

A simple pole at the origin

```mathematica
In[10]:= Residue[1/z, {z, 0}]
Out[10]= 1
```

A double pole has no (z)^-1 term, so residue 0

```mathematica
In[11]:= Residue[1/z^2, {z, 0}]
Out[11]= 0
```

A transcendental integrand expanded directly

```mathematica
In[12]:= Residue[Cot[z], {z, 0}]
Out[12]= 1
```

A simple pole at a complex location

```mathematica
In[13]:= Residue[1/(z^2 + 1), {z, I}]
Out[13]= -1/2*I
```

Residue 1 at the simple pole of Zeta

```mathematica
In[14]:= Residue[Zeta[s], {s, 1}]
Out[14]= 1
```

A branch point -- undefined, left unevaluated

```mathematica
In[15]:= Residue[1/Sqrt[z], {z, 0}]
Out[15]= Residue[1/Sqrt[z], {z, 0}]
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

**Algorithm.** `builtin_residue` handles `Residue[f, {z, z0}]` by reading the
coefficient of `(z - z0)^-1` out of the Laurent expansion of `f`.
`residue_compute` first classifies the integrand with `residue_is_rational_in`
(is `Together[f]` a ratio of polynomials in `z`?). For a rational `f` it tries a
simple-pole fast path, `residue_simple_pole`: with `P/Q = Together[f]` and `Q`
having a simple zero at `z0` (`Q(z0) == 0` via `PossibleZeroQ`, `Q'(z0) != 0`),
the residue is `P(z0)/Q'(z0)`, which bypasses the series inverter entirely.
Failing that, `residue_shift_form` substitutes `z -> z0 + w` and `Expand`s the
shifted numerator and denominator *separately* (keeping them coprime), then
`residue_extract` runs `Series[..., {w, 0, 0}]` and reads the `w^-1` term.
Transcendental integrands (`Cot`, `Zeta` near its pole, unknown `f[z]`) skip the
shift and expand directly about `z0`, so the series engine can use its built-in
knowledge of the function's Laurent series there.

**Data structures.** `Expr` trees throughout, driven by `eval_and_free` wrappers
(`residue_eval1`/`residue_eval2`) that build and evaluate `Together`,
`Numerator`, `Denominator`, `Expand`, `D`, `ReplaceAll`, and `Series`. The
expansion returns a `SeriesData[z, z0, {coefs}, nmin, nmax, den]`;
`residue_extract` reads the coefficient at index `-1 - nmin`. A fresh local
`w = Residue\`$w` carries the shifted expansion.

**Complexity / limits.** The simple-pole path is one differentiation plus two
substitutions; the series path costs a Laurent expansion whose order
`residue_extract` raises adaptively (capped at 256) until the `-1` term is
resolved. A residue is defined only for an ordinary Laurent expansion
(`den == 1`): a fractional-power (Puiseux) expansion signals a branch point, so
`Residue[1/Sqrt[z], {z, 0}]` is left unevaluated. Only the two-argument form is
handled; fewer than two arguments emit `Residue::argm`. The separate-shift
design is what keeps an algebraic pole location (a `z0` that is a sum of
radicals) tractable, since re-cancelling the shifted ratio can spawn a second
spelling of the same algebraic number and blow the series inversion up
combinatorially.

**Attributes:** `Protected`.

## References

**See also:** [NResidue](../../numerical-calculus/NResidue/), [Together](../../algebra/Together/), [Zeta](../../special-functions/Zeta/)

- Source: [`src/calculus/residue.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/residue.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_residue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_residue.c)

## Notes & additional examples

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
