# $MaxNumber

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$MaxNumber
    gives the maximum arbitrary-precision number that can be
    represented on this computer system.

With USE_MPFR builds, this is the largest finite value at machine
precision under MPFR's current exponent range; otherwise it equals
$MaxMachineNumber.
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

A Protected OwnValue registered in `system_constants_init` (`src/core.c`). In a `USE_MPFR` build it is the largest finite value at machine precision (`DBL_MANT_DIG` bits), computed by `mpfr_set_inf` then `mpfr_nextbelow` and stored via `expr_new_mpfr_move`; without MPFR it collapses to `expr_new_real(DBL_MAX)`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
