### Worked examples

```mathematica
In[1]:= Head[$MaxNumber]  (* a Real; its literal value is build-dependent *)
```

```mathematica
In[1]:= $MaxNumber >= $MaxMachineNumber  (* never smaller than the machine maximum *)
```

### Notes

`$MaxNumber` is a read-only Protected OwnValue: the largest finite number Mathilda can
represent. In a `USE_MPFR` build it is the largest finite value at machine precision
(`DBL_MANT_DIG` bits) with MPFR's enormous exponent range — computed by `mpfr_set_inf`
then `mpfr_nextbelow` — so it is astronomically larger than `$MaxMachineNumber`. Without
MPFR there is no arbitrary-precision representation and it collapses onto `DBL_MAX`,
equal to `$MaxMachineNumber`.

Because the literal value therefore depends on whether MPFR was linked, the examples above
are structural; `$MaxNumber >= $MaxMachineNumber` holds in either build.
