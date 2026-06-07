# $MinMachineNumber

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$MinMachineNumber
    gives the smallest positive machine-precision number that can
    be represented in normalized form on this computer system.

Equals the platform's DBL_MIN.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= $MachinePrecision
Out[1]= 15.9546

In[2]:= $MachineEpsilon
Out[2]= 2.22045e-16

In[3]:= {$MinMachineNumber, $MaxMachineNumber}
Out[3]= {2.22507e-308, 1.79769e+308}

In[4]:= MachineNumberQ[$MaxNumber]   (* MPFR, not machine *)
Out[4]= False
```

## Implementation notes

A Protected OwnValue set in `system_constants_init` (`src/core.c`) to `expr_new_real(DBL_MIN)`, the smallest positive normalized IEEE 754 `double`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
