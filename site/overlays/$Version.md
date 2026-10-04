### Worked examples

```mathematica
In[1]:= StringQ[$Version]  (* a descriptive string, not a number *)
```

```mathematica
In[1]:= Head[$Version]  (* its head is String *)
```

### Notes

`$Version` is a read-only Protected OwnValue holding a descriptive string, e.g.
`"Mathilda 0.266 (<compiler>, GMP ..., MPFR ..., FLINT ..., ...)"`, assembled at compile
time from `MATHILDA_VERSION_STRING` and each linked library's own version macros — so it
names exactly what this binary was built against. The numeric release is `$VersionNumber`.

The exact text changes every release, so the examples check it structurally rather than
against a literal.
