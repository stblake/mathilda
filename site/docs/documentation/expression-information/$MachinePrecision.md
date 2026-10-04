# $MachinePrecision

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$MachinePrecision`**

gives the number of decimal digits of precision used for machine-precision numbers.

<details>
<summary>Notes</summary>

Derived from the platform's DBL\_MANT\_DIG -- typically 53\*Log\[10,2\] (~ 15.9546) on IEEE 754 systems.

</details>

## Examples (5)

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

### Applications (1)

Decimal digits carried by a machine double, 53 log10 2

```mathematica
In[5]:= $MachinePrecision
Out[5]= 15.9546
```

## Implementation notes

A Protected OwnValue registered in `system_constants_init` (`src/core.c`) as `expr_new_real(NUMERIC_MACHINE_PRECISION_DIGITS)` — the number of decimal digits in a machine `double` (~15.95).

**Attributes:** `Protected`.

## References

**See also:** [$MachineEpsilon](../../expression-information/$MachineEpsilon/), [$MinMachineNumber](../../expression-information/$MinMachineNumber/), [$MaxMachineNumber](../../expression-information/$MaxMachineNumber/), [$MaxNumber](../../expression-information/$MaxNumber/), [$MinNumber](../../expression-information/$MinNumber/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Notes

`$MachinePrecision` is a read-only Protected OwnValue giving the number of decimal digits
in a machine `double`, `53 · log₁₀ 2 ≈ 15.9546`. It is the precision `N[expr]` (with no
second argument) works at, and the threshold that separates machine-precision reals from
the arbitrary-precision MPFR path.

The value is a fixed machine constant on every IEEE 754 platform.
