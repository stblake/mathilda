### Worked examples

```mathematica
(* HoldAll and Protected: the body is held until each frame is drawn *)
In[1]:= Attributes[Animate]
```

```mathematica
(* a rotating sine wave; the window opens interactively and blocks the REPL *)
In[1]:= Animate[Plot[Sin[x + t], {x, 0, 2 Pi}], {t, 0, 2 Pi}]
```

```mathematica
(* a bouncing disk, ping-ponging over the range *)
In[1]:= Animate[Graphics[Disk[{t, 0}, 0.5], PlotRange -> {{0, 5}, {-1, 1}}], {t, 0, 5}, AnimationDirection -> ForwardBackward]
```

```mathematica
(* 3D content gets an orbit camera on top of the playback controls *)
In[1]:= Animate[Plot3D[Sin[x + t] Cos[y], {x, -3, 3}, {y, -3, 3}], {t, 0, 2 Pi}]
```

### Notes

`Animate` is `HoldAll` and `Protected`: the body `expr` is held and re-evaluated
every frame with the iterator variable bound to the current parameter value, a
single phase `phi` driving all iterators. It opens a Raylib window with playback
controls — `AnimationDirection`, `AnimationRate`, `AnimationRepetitions`,
`AnimationRunning`, `DefaultDuration`, `ControlPlacement`, `RefreshRate` — and the
event loop **blocks the REPL** until the window closes, whereupon it returns `Null`.
A `Graphics3D`/`Plot3D` body gets its own persistent orbit camera.

Under headless conditions (the `MATHILDA_NO_WINDOW` environment variable, as used by
documentation generation) `Animate` returns unevaluated — the body is never run and
no window opens — which is why the examples above display as the held
`Animate[...]` form here rather than opening a window.
