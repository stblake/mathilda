# LogIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LogIntegral[z]
    gives the logarithmic integral li(z), the principal value of
    Integral_0^z dt/ln t, equal to ExpIntegralEi[Log[z]], with a branch cut
    on (-Infinity, 1). LogIntegral[0] = 0, LogIntegral[1] = -Infinity,
LogIntegral[Infinity] = Infinity. Real and complex inputs evaluate
numerically at machine or arbitrary (MPFR) precision; D[LogIntegral[z], z] =
1/Log[z]. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact special values: `LogIntegral[0] = 0`, `LogIntegral[1] = -Infinity`,

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= N[LogIntegral[2], 40]
Out[1]= 1.0451637801174927848445888891946131365227
```

```mathematica
In[1]:= D[LogIntegral[z], z]
Out[1]= 1/Log[z]
```

```mathematica
In[1]:= N[LogIntegral[10^6], 30]
Out[1]= 78627.54915946218191986291074769
```

```mathematica
In[1]:= N[LogIntegral[1000], 20]
Out[1]= 177.609657990152226688
```

### Notes

`LogIntegral[z]` is the logarithmic integral `li(z)`, the principal value of `Integral_0^z dt/Log[t]`, equal to `ExpIntegralEi[Log[z]]`, with a branch cut on `(-Infinity, 1)`. Its derivative is `1/Log[z]`. `li(x)` is the leading term of the prime-counting approximation `PrimePi[x] ~ li(x)`; for example `li(10^6)` is about `78627.5`, close to `PrimePi[10^6] = 78498`. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision. `LogIntegral[1] = -Infinity`. Listable.
