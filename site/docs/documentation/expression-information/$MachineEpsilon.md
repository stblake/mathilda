# $MachineEpsilon

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$MachineEpsilon`**

gives the difference between 1.0 and the next-nearest number representable as a machine-precision number.

<details>
<summary>Notes</summary>

Equals the platform's DBL\_EPSILON; measures the granularity of machine-precision numbers.

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

The smallest e for which 1 + e differs from 1 in machine arithmetic

```mathematica
In[5]:= $MachineEpsilon
Out[5]= 2.22045e-16
```

Exactly two to the minus fifty-two for an IEEE 754 double

```mathematica
In[6]:= $MachineEpsilon == 2^-52
Out[6]= True
```

## Implementation notes

A read-only system constant bound as an OwnValue in `system_constants_init` (`src/core.c`) via `register_system_constant`, then marked `ATTR_PROTECTED`. Its value is `expr_new_real(DBL_EPSILON)` — the `<float.h>` machine epsilon of the local IEEE 754 `double`.

**Attributes:** `Protected`.

## References

**See also:** [$MachinePrecision](../../expression-information/$MachinePrecision/), [$MinMachineNumber](../../expression-information/$MinMachineNumber/), [$MaxMachineNumber](../../expression-information/$MaxMachineNumber/), [$MaxNumber](../../expression-information/$MaxNumber/), [$MinNumber](../../expression-information/$MinNumber/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Notes

`$MachineEpsilon` is a read-only system constant: the `<float.h>` `DBL_EPSILON` of the
local IEEE 754 `double`, bound as a Protected OwnValue. It is the gap between `1.0` and
the next representable double, `2^-52 ≈ 2.22045×10⁻¹⁶`, and is the natural unit for
machine-precision round-off tolerances.

The value is a fixed machine constant on every IEEE 754 platform, so it is safe to quote
literally.
