### Worked examples

```mathematica
(* HoldAll and Protected: the body is held until a control is drawn *)
In[1]:= Attributes[Manipulate]
```

```mathematica
(* a continuous slider driving the frequency *)
In[1]:= Manipulate[Plot[Sin[n x], {x, 0, 2 Pi}], {n, 1, 5}]
```

```mathematica
(* a discrete button set switches between functions *)
In[1]:= Manipulate[Plot[f, {x, -5, 5}], {f, {Sin[x], Cos[x], x^2}}]
```

```mathematica
(* explicit default and step: {{var, default}, min, max, step} *)
In[1]:= Manipulate[Graphics[Disk[{0, 0}, r], PlotRange -> {{-5, 5}, {-5, 5}}], {{r, 2}, 0.5, 5, 0.25}]
```

### Notes

`Manipulate` is `HoldAll` and `Protected`: one control row is drawn per variable and
the body is re-evaluated with each variable bound to its current value as a control
changes. Control specs cover continuous sliders (`{u, min, max}`, with an optional
step `du`, or a `{{u, u0}, ...}` explicit default) and discrete button sets
(`{u, {v1, v2, ...}}`). A footer **Reset** button restores every control's default.
Unlike `Animate` there is no playback transport — every control is user-driven.

The window blocks the REPL until closed, then returns `Null`; a `Graphics3D`/`Plot3D`
body gets its own orbit camera. Under the headless `MATHILDA_NO_WINDOW` environment
(used by documentation generation) `Manipulate` returns unevaluated, which is why the
examples above display as the held `Manipulate[...]` form rather than opening a
window.
