# $MaxMachineNumber

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$MaxMachineNumber`**

gives the largest machine-precision number that can be used on this computer system.

<details>
<summary>Notes</summary>

Equals the platform's DBL\_MAX.

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
```

MPFR, not machine

```mathematica
In[4]:= MachineNumberQ[$MaxNumber]
Out[4]= False
```

## Implementation notes

A Protected OwnValue set in `system_constants_init` (`src/core.c`) to `expr_new_real(DBL_MAX)`, the largest finite IEEE 754 `double`.

**Attributes:** `Protected`.

## References

**See also:** [$MachinePrecision](../../expression-information/$MachinePrecision/), [$MachineEpsilon](../../expression-information/$MachineEpsilon/), [$MinMachineNumber](../../expression-information/$MinMachineNumber/), [$MaxNumber](../../expression-information/$MaxNumber/), [$MinNumber](../../expression-information/$MinNumber/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
