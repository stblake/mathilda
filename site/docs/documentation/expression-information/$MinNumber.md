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

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= $MachinePrecision
Out[1]= 15.9546

In[2]:= $MachineEpsilon
Out[2]= 2.22045e-16

In[3]:= {$MinMachineNumber, $MaxMachineNumber}
Out[3]= {2.22507e-308, 1.79769e+308}
```

MPFR, not machine

```mathematica
In[4]:= MachineNumberQ[$MaxNumber]
Out[4]= False
```

### Applications (2)

A Real; its literal value is build-dependent

```mathematica
In[5]:= Head[$MinNumber]
Out[5]= Real
```

Never larger than the machine minimum

```mathematica
In[6]:= $MinNumber <= $MinMachineNumber
Out[6]= True
```

## Implementation notes

A Protected OwnValue registered in `system_constants_init` (`src/core.c`). Under `USE_MPFR` it is the smallest positive value at machine precision (`mpfr_set_zero` then `mpfr_nextabove`, stored via `expr_new_mpfr_move`); without MPFR it collapses to `expr_new_real(DBL_MIN)`.

**Attributes:** `Protected`.

## References

**See also:** [$MachinePrecision](../../expression-information/$MachinePrecision/), [$MachineEpsilon](../../expression-information/$MachineEpsilon/), [$MinMachineNumber](../../expression-information/$MinMachineNumber/), [$MaxMachineNumber](../../expression-information/$MaxMachineNumber/), [$MaxNumber](../../expression-information/$MaxNumber/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Notes

`$MinNumber` is a read-only Protected OwnValue: the smallest positive number Mathilda can
represent. In a `USE_MPFR` build it is the smallest positive value at machine precision
(`mpfr_set_zero` then `mpfr_nextabove`), far below `$MinMachineNumber` thanks to MPFR's
exponent range. Without MPFR it collapses onto `DBL_MIN`, equal to `$MinMachineNumber`.

Because the literal value therefore depends on whether MPFR was linked, the examples above
are structural; `$MinNumber <= $MinMachineNumber` holds in either build.
