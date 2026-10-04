### Worked examples

```mathematica
In[1]:= $MachineEpsilon  (* the smallest e for which 1 + e differs from 1 in machine arithmetic *)
```

```mathematica
In[1]:= $MachineEpsilon == 2^-52  (* exactly two to the minus fifty-two for an IEEE 754 double *)
```

### Notes

`$MachineEpsilon` is a read-only system constant: the `<float.h>` `DBL_EPSILON` of the
local IEEE 754 `double`, bound as a Protected OwnValue. It is the gap between `1.0` and
the next representable double, `2^-52 ≈ 2.22045×10⁻¹⁶`, and is the natural unit for
machine-precision round-off tolerances.

The value is a fixed machine constant on every IEEE 754 platform, so it is safe to quote
literally.
