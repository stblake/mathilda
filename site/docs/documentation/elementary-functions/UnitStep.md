# UnitStep

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
UnitStep[x]
    gives 0 for x < 0 and 1 for x >= 0 (the value at 0 is 1).
UnitStep[x1, x2, ...]
    gives 1 only when none of the xi are negative, otherwise 0.
UnitStep[] is 1. The result is always exact. Exact symbolic real
arguments are resolved by numerical certification; non-real or
unresolved arguments are left unevaluated. UnitStep is Listable and
Orderless.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= UnitStep[0]
Out[1]= 1

In[2]:= UnitStep[1, Pi, 5.3]
Out[2]= 1

In[3]:= UnitStep[{-1.6, 3.200000000000}]
Out[3]= {0, 1}

In[4]:= UnitStep[Sqrt[2] - 99/70]
Out[4]= 0

In[5]:= D[UnitStep[x], x]
Out[5]= Piecewise[{{Indeterminate, x == 0}}, 0]

In[6]:= D[UnitStep[x, y, z], z]
Out[6]= UnitStep[x, y] Piecewise[{{Indeterminate, z == 0}}, 0]
```

## Implementation notes

- `Listable`, `NumericFunction`, `Orderless`, `Protected`.
- The result is **always exact** -- an integer `0` or `1` -- for real numeric

**Attributes:** `Listable`, `NumericFunction`, `Orderless`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= UnitStep[-2]
Out[1]= 0

In[2]:= UnitStep[3]
Out[2]= 1
```

UnitStep is `Listable`, so it vectorizes over a table to produce a discrete step profile (here switching on at `k = 3`):

```mathematica
In[1]:= Table[UnitStep[k - 3], {k, 0, 6}]
Out[1]= {0, 0, 0, 1, 1, 1, 1}
```

Exact symbolic-real arguments are resolved by numerical certification — even transcendental comparisons collapse to an exact `0` or `1`:

```mathematica
In[1]:= UnitStep[Pi - 3]
Out[1]= 1

In[2]:= UnitStep[Log[2] - Log[3]]
Out[2]= 0
```

The multivariate form is the indicator of the nonnegative orthant, returning `1` only when no argument is negative:

```mathematica
In[1]:= UnitStep[1, -1, 2]
Out[1]= 0
```

### Notes

`UnitStep[x]` is `0` for `x < 0` and `1` for `x >= 0` (the value at `0` is `1`). The result is always exact: certifiable real arguments resolve numerically, while non-real or unresolved arguments are left unevaluated. UnitStep is `Listable` and `Orderless`.
