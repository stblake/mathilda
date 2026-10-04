### Worked examples

```mathematica
In[1]:= Head[$VersionNumber]  (* its head is Real *)
```

```mathematica
In[1]:= $VersionNumber > 0  (* a positive release number *)
```

### Notes

`$VersionNumber` is a read-only Protected OwnValue giving the Mathilda release as a Real,
such as `0.266`. It comes straight from `MATHILDA_VERSION_NUMBER` in `src/version.h`, the
single source of truth for the release; the descriptive form is `$Version`.

The numeric value changes every release (and the printer drops trailing zeros, so `0.160`
would print as `0.16`), so the examples test it structurally rather than against a literal.
