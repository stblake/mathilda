# $MinNumber

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$MinNumber`**

gives the minimum positive arbitrary-precision number that can be represented on this computer system.

<details>
<summary>Notes</summary>

With USE\_MPFR builds, this is the smallest positive value at machine precision under MPFR's current exponent range; otherwise it equals $MinMachineNumber.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

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

A Protected OwnValue registered in `system_constants_init` (`src/core.c`). Under `USE_MPFR` it is the smallest positive value at machine precision (`mpfr_set_zero` then `mpfr_nextabove`, stored via `expr_new_mpfr_move`); without MPFR it collapses to `expr_new_real(DBL_MIN)`.

**Attributes:** `Protected`.

## See also

[$MachinePrecision](../../expression-information/$MachinePrecision/), [$MachineEpsilon](../../expression-information/$MachineEpsilon/), [$MinMachineNumber](../../expression-information/$MinMachineNumber/), [$MaxMachineNumber](../../expression-information/$MaxMachineNumber/), [$MaxNumber](../../expression-information/$MaxNumber/)

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
