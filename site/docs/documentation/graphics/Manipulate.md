# Manipulate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Manipulate[expr, {u, umin, umax}, ...]`**

Opens an interactive window with one control per variable and re-evaluates expr with each variable bound to its current value whenever a control changes. Returns Null once the window is closed. expr is typically a Graphics\[...\] or Plot\[...\] call that depends on the control variables. Unlike Animate, every control is independently user-driven from the first frame -- there is no animation phase or playback transport. Control specs (any number, one row each): {u, umin, umax}              continuous slider, default = umin {u, umin, umax, du}          continuous slider with step du {{u, u0}, umin, umax}        continuous slider, explicit default u0 {{u, u0}, umin, umax, du}    continuous slider, default + step {u, {v1, v2, ...}}           discrete button set, default = v1 {{u, u0}, {v1, v2, ...}}     discrete button set, explicit default u0 Click-drag a slider handle or click a discrete button to change its value; click Reset to restore every control to its default. Esc closes the window.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Manipulate[Plot[Sin[n x], {x, 0, 2 Pi}], {n, 1, 5}]
Out[1]= Manipulate[Plot[Sin[n x], {x, 0, 2 Pi}], {n, 1, 5}]

In[2]:= Manipulate[Plot[f, {x, -5, 5}], {f, {Sin[x], Cos[x], x^2}}]
Out[2]= Manipulate[Plot[f, {x, -5, 5}], {f, {Sin[x], Cos[x], x^2}}]

In[3]:= Manipulate[Plot[a Sin[x] + b, {x, 0, 2 Pi}], {a, 0, 3}, {b, {-1, 0, 1}}]
Out[3]= Manipulate[Plot[a Sin[x] + b, {x, 0, 2 Pi}], {a, 0, 3}, {b, {-1, 0, 1}}]

In[4]:= Manipulate[Plot3D[Sin[x + n] Cos[y], {x, -3, 3}, {y, -3, 3}], {n, 0, 3}]
Out[4]= Manipulate[Plot3D[Sin[x + n] Cos[y], {x, -3, 3}, {y, -3, 3}], {n, 0, 3}]
```

### Options (1)

```mathematica
In[5]:= Manipulate[ Graphics[Disk[{0, 0}, r], PlotRange -> {{-5, 5}, {-5, 5}}], {{r, 2}, 0.5, 5, 0.25}]
Out[5]= Manipulate[-Graphics-, {{r, 2}, 0.5, 5, 0.25}]
```

## Algorithm

manipulate.c — Manipulate[expr, {u, umin, umax}, ...]

Opens a Raylib window (when USE_GRAPHICS is compiled in) with one control row per variable, stacked above a thin footer:

```text
  ┌─────────────────────────────────────────────────────────────┐
  │                                                               │
  │                        (rendered content)                    │
  │                                                               │
  ├─────────────────────────────────────────────────────────────┤
  │ u  │ 0 ─────────────●────────── 10              │  3.42     │
  │ f  │ [Sin[x]] [Cos[x]] [Tan[x]]                              │
  ├─────────────────────────────────────────────────────────────┤
  │ [Reset]  Esc: close                                          │
  └─────────────────────────────────────────────────────────────┘
```

Each control is independently user-driven from the first frame — there is no animation phase or playback transport (see Animate for that). A control is either:

```text
  - a continuous drag slider:  {u, umin, umax}, {u, umin, umax, du},
    {{u, u0}, umin, umax}, {{u, u0}, umin, umax, du}
  - a discrete button set:     {u, {v1, v2, ...}}, {{u, u0}, {v1, ...}}
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Animate](../../graphics/Animate/), [Show](../../graphics/Show/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_manipulate.c`](https://github.com/stblake/mathilda/blob/main/tests/test_manipulate.c)
