### Worked examples

```mathematica
In[1]:= Green  (* the named primary evaluates straight to its RGBColor literal *)
```

```mathematica
In[1]:= FullForm[Green]  (* an OwnValue, not an inert head -- it has already resolved *)
```

```mathematica
In[1]:= Head[Graphics[{Green, Line[{{0, 0}, {1, 1}}]}]]  (* usable wherever a directive is *)
```

### Notes

`Green` is the named colour `RGBColor[0, 1, 0]`. It is bound as a `Protected` OwnValue, so
it resolves to a real colour literal wherever a bare argument is evaluated: `FullForm[Green]`
is `RGBColor[0, 1, 0]` and `Head[Green]` is `RGBColor`.

That resolution is what lets it be used anywhere a style directive is accepted — inside
`Graphics`, as a `PlotStyle`/`Epilog` value, or in a `Directive`. The exact-`0`/`1`
components print as integers. `Green` carries no opacity; use `Opacity[a]` or the
four-argument `RGBColor[r, g, b, a]`.
