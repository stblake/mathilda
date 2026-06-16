# PolyLog

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PolyLog[n, z]
    gives the polylogarithm Li_n(z) = Sum_{k>=1} z^k/k^n.
PolyLog[n, p, z]
    gives the Nielsen generalized polylogarithm S_{n,p}(z) (accepted but
left unevaluated).
Special arguments reduce in closed form: PolyLog[n, 0] = 0,
PolyLog[1, z] = -Log[1-z], PolyLog[0, z] = z/(1-z), negative integer
orders give rational functions, PolyLog[n, 1] = Zeta[n] and
PolyLog[n, -1] = (2^(1-n)-1) Zeta[n] for integer n >= 2, with exact forms
for PolyLog[2, 1/2] and PolyLog[3, 1/2]. Inexact real or complex arguments
evaluate numerically at machine or arbitrary (MPFR) precision via a power
series or the Jonquiere/zeta expansion. There is a branch cut from 1 to
Infinity. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PolyLog[3, 1/2]
Out[1]= 1/6 Log[2]^3 - 1/12 Log[2] Pi^2 + 7/8 Zeta[3]

In[2]:= PolyLog[2, 0.9]
Out[2]= 1.29971
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

At `z = 1` the polylogarithm is the Riemann zeta value; `PolyLog[1, z]` is the
ordinary logarithm:

```mathematica
In[1]:= PolyLog[2, 1]
Out[1]= 1/6 Pi^2

In[2]:= PolyLog[1, z]
Out[1]= -Log[1 - z]
```

The classical closed forms at `z = 1/2` — including the celebrated Euler value
of `Li_3(1/2)` mixing `Log[2]`, `Pi^2`, and `Zeta[3]`:

```mathematica
In[1]:= PolyLog[2, 1/2]
Out[1]= -1/2 Log[2]^2 + 1/12 Pi^2

In[2]:= PolyLog[3, 1/2]
Out[2]= 1/6 Log[2]^3 - 1/12 Log[2] Pi^2 + 7/8 Zeta[3]
```

Non-positive integer orders are rational functions of `z` (the negative-order
polylogarithms), produced in closed form:

```mathematica
In[1]:= PolyLog[0, z]
Out[1]= z/(1 - z)

In[2]:= PolyLog[-2, z]
Out[2]= (z + z^2)/(1 - z)^3
```

The numeric core evaluates real and complex arguments to arbitrary precision —
here `Li_2(1/2)` to 40 digits and `Li_3` at a complex point to 30:

```mathematica
In[1]:= N[PolyLog[2, 1/2], 40]
Out[1]= 0.58224052646501250590265632015968010874412

In[2]:= N[PolyLog[3, 1/2 + I/2], 30]
Out[2]= 0.48615953708556007896672148708 + 0.5700774070887689781956097575898*I
```

### Notes

`PolyLog[n, z]` is the polylogarithm Li_n(z) = Σ_{k≥1} z^k / k^n. Special
arguments reduce in closed form: `PolyLog[n, 0] = 0`, `PolyLog[1, z] = -Log[1-z]`,
`PolyLog[0, z] = z/(1-z)`, negative integer orders give rational functions, and
`PolyLog[n, 1] = Zeta[n]`, `PolyLog[n, -1] = (2^(1-n)-1) Zeta[n]` for integer
`n ≥ 2`, with exact forms for `PolyLog[2, 1/2]` and `PolyLog[3, 1/2]`. Inexact
real or complex arguments evaluate numerically at machine or MPFR precision via
a power series or the Jonquière/zeta expansion. There is a branch cut from `1`
to `Infinity`. Listable. The Nielsen generalized form `PolyLog[n, p, z]` is
accepted but left unevaluated.
