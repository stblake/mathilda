# $MachinePrecision

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$MachinePrecision
    gives the number of decimal digits of precision used for
    machine-precision numbers.

Derived from the platform's DBL_MANT_DIG -- typically 53*Log[10,2]
(~ 15.9546) on IEEE 754 systems.
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

A Protected OwnValue registered in `system_constants_init` (`src/core.c`) as `expr_new_real(NUMERIC_MACHINE_PRECISION_DIGITS)` — the number of decimal digits in a machine `double` (~15.95).

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
