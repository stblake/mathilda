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
