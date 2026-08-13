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

A read-only system constant bound as an OwnValue in `system_constants_init` (`src/core.c`) via `register_system_constant`, then marked `ATTR_PROTECTED`. Its value is `expr_new_real(DBL_EPSILON)` — the `<float.h>` machine epsilon of the local IEEE 754 `double`.

**Attributes:** `Protected`.

## References

**See also:** [$MachinePrecision](../../expression-information/$MachinePrecision/), [$MinMachineNumber](../../expression-information/$MinMachineNumber/), [$MaxMachineNumber](../../expression-information/$MaxMachineNumber/), [$MaxNumber](../../expression-information/$MaxNumber/), [$MinNumber](../../expression-information/$MinNumber/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
