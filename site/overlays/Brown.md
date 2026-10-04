### Worked examples

```mathematica
In[1]:= Brown  (* a tertiary tone -- components are reals, not 0/1 *)
```

```mathematica
In[1]:= FullForm[Brown]  (* an OwnValue, not an inert head -- already resolved *)
```

```mathematica
In[1]:= Head[Graphics[{Brown, Disk[]}]]  (* usable wherever a directive is expected *)
```

### Notes

`Brown` is the named colour `RGBColor[0.6, 0.4, 0.2]`. It is bound as a `Protected`
OwnValue, so it resolves to a real colour literal wherever a bare argument is evaluated:
`FullForm[Brown]` is `RGBColor[0.6, 0.4, 0.2]` and `Head[Brown]` is `RGBColor`.

Its components are machine reals (none is exactly `0` or `1`), so they print as reals,
unlike the pure primaries. The light companion `LightBrown` (`RGBColor[0.94, 0.91, 0.88]`)
is a separate constant, not a programmatic lightening of `Brown`. `Brown` carries no
opacity; use `Opacity[a]` or the four-argument `RGBColor[r, g, b, a]`.
