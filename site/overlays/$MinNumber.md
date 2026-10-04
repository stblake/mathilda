### Worked examples

```mathematica
In[1]:= Head[$MinNumber]  (* a Real; its literal value is build-dependent *)
```

```mathematica
In[1]:= $MinNumber <= $MinMachineNumber  (* never larger than the machine minimum *)
```

### Notes

`$MinNumber` is a read-only Protected OwnValue: the smallest positive number Mathilda can
represent. In a `USE_MPFR` build it is the smallest positive value at machine precision
(`mpfr_set_zero` then `mpfr_nextabove`), far below `$MinMachineNumber` thanks to MPFR's
exponent range. Without MPFR it collapses onto `DBL_MIN`, equal to `$MinMachineNumber`.

Because the literal value therefore depends on whether MPFR was linked, the examples above
are structural; `$MinNumber <= $MinMachineNumber` holds in either build.
