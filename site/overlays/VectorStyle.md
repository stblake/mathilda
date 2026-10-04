### Worked examples

```mathematica
In[1]:= Head[VectorStyle]  (* a bare VectorPlot option symbol *)
```

```mathematica
In[1]:= Attributes[VectorStyle]  (* inert and Protected -- no builtin behaviour *)
```

```mathematica
In[1]:= FullForm[VectorStyle -> Red]  (* the named colour resolves to RGBColor[1, 0, 0] *)
```

### Notes

`VectorStyle` is an option for `VectorPlot` that applies style directive(s) globally to
every arrow. It is an inert, `Protected` option-name symbol — no builtin, no DownValues —
so it does nothing on its own; the `VectorPlot` renderer reads it from the call's option
list.

Its value is one or more directives (`RGBColor`, `Thickness`, …); a named colour resolves,
so `VectorStyle -> Red` is stored as `Rule[VectorStyle, RGBColor[1, 0, 0]]`. It overrides
any per-arrow `ColorFunction` and affects only the arrows' appearance, never the sampled
field.
