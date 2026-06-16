# LogGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LogGamma[z]
    gives the log-gamma function log(Gamma(z)), analytic except for a branch
cut on the negative reals. Exact at integer and half-integer z (with the
negative-axis branch term), divergent (Infinity) at non-positive integers,
and evaluated numerically for real and complex z at machine or arbitrary
(MPFR) precision. D[LogGamma[z], z] is PolyGamma[0, z]. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- **Exact closed forms.** Integers reduce as `LogGamma[n] = Log[(n-1)!]`

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= LogGamma[5]
Out[1]= Log[24]
```

```mathematica
In[1]:= LogGamma[1/2]
Out[1]= Log[Sqrt[Pi]]
```

```mathematica
In[1]:= D[LogGamma[z], z]
Out[1]= PolyGamma[0, z]
```

```mathematica
In[1]:= N[LogGamma[100], 40]
Out[1]= 359.13420536957539877604401046028690961264
```

```mathematica
In[1]:= N[LogGamma[1 + I], 30]
Out[1]= -0.6509231993018563388852168315042 - 0.3016403204675331978875316577968*I
```

### Notes

`LogGamma[z]` is `log(Gamma(z))`, analytic except for a branch cut on the negative reals. It is exact at integer and half-integer arguments (`LogGamma[5]` is `Log[4!] = Log[24]`), divergent at non-positive integers, and evaluates numerically for real or complex `z` at machine or arbitrary (MPFR) precision. Its derivative is `PolyGamma[0, z]`. Unlike `Log[Gamma[z]]`, `LogGamma` tracks the correct sheet, which matters for large or complex arguments where `Gamma` overflows. Listable.
