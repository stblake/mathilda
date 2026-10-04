### Worked examples

```mathematica
(* accepted on PolarPlot; Cartesian axes are drawn until the polar grid lands *)
In[1]:= PolarPlot[2, {t, 0, 2 Pi}, PolarAxes -> True]
```

```mathematica
(* the option rides through verbatim as a Rule on the Graphics[] object *)
In[1]:= Cases[PolarPlot[Sin[2 t], {t, 0, 2 Pi}, PolarAxes -> True], (PolarAxes -> v_) :> v, Infinity]
```

```mathematica
(* an inert, Protected option keyword *)
In[1]:= MemberQ[Attributes[PolarAxes], Protected]
```

```mathematica
(* it has no values of its own, so it stays symbolic *)
In[1]:= Head[PolarAxes]
```

### Notes

`PolarAxes` is not a plotter but an inert, `Protected` **option keyword** for
`PolarPlot`. `PolarAxes -> True` asks for a polar grid overlay — radial circles at
regular intervals plus angular degree/radian labels — but that overlay is not yet
rendered, so the plot is drawn with ordinary Cartesian axes instead. This is a
documented placeholder rather than a silent no-op: the option is recognised and
preserved (it passes through as a `Rule` onto the returned `Graphics[...]`), so
existing `PolarPlot` calls that set it keep working unchanged.

The symbol carries no builtin and no own-values; it has only the `Protected`
attribute and a docstring, and stays symbolic in every context.
