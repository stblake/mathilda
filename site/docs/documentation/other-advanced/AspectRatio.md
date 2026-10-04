# AspectRatio

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AspectRatio`**

is an option for Graphics and Plot that specifies the ratio of height to width of the rendered plot.

<details>
<summary>Notes</summary>

AspectRatio -\> Automatic sets the ratio from the actual coordinate values (true geometry); AspectRatio -\> Full stretches the graphics to fill the enclosing region; AspectRatio -\> a uses the explicit height-to-width ratio a. Plot defaults to 1/GoldenRatio.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

A bare option symbol, not a function

```mathematica
In[1]:= Head[AspectRatio]
Out[1]= Symbol
```

It only ever sits on the left of a rule

```mathematica
In[2]:= FullForm[AspectRatio -> 1]
Out[2]= Rule[AspectRatio, 1]
```

Plot's stored default is Automatic

```mathematica
In[3]:= Options[Plot, AspectRatio]
Out[3]= {AspectRatio -> Automatic}
```

Accepted, still a Graphics

```mathematica
In[4]:= Head[Graphics[{Disk[]}, AspectRatio -> 1, ImageSize -> 300]]
Out[4]= Graphics
```

## Implementation notes

**Definition.** `AspectRatio` is an **option name** for `Graphics` and `Plot` giving the
ratio of rendered height to width. It is a bare option symbol, not a function: it has no
builtin, no DownValues and (unusually) not even the `Protected` attribute — its only
registration is the docstring set in `info_init` (`src/info.c`). Its meaning lives entirely
in the graphics back end (`src/graphics/`), which reads the option out of the option list
a `Graphics`/`Plot` call carries.

**Representation.** `AspectRatio` stays an inert `EXPR_SYMBOL`; it appears only on the
left of a rule, so `FullForm[AspectRatio -> 1]` is `Rule[AspectRatio, 1]`. The value may
be an explicit height-to-width real, `Automatic` (derive the ratio from the actual
coordinate values — true geometry), or `Full` (stretch to fill the enclosing region).
`Options[Plot, AspectRatio]` reports `Plot`'s stored default, `Automatic`; a standalone
`Plot` curve otherwise lays out at `1/GoldenRatio`.

**Usage & limits.** Supplied as `Graphics[prims, AspectRatio -> r]` or
`Plot[f, {x, a, b}, AspectRatio -> r]`, it is consumed by the renderer when the window or
export is laid out; it changes the displayed shape, never the plotted data, so a call that
accepts it still reduces to an ordinary `Graphics[...]` object. Because the symbol is not
`Protected`, it could in principle be assigned a value — it is intended only as an option
tag.

**Attributes:** none registered.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`AspectRatio` is an option for `Graphics` and `Plot` that fixes the ratio of rendered
height to width. It is an inert option-name symbol — no builtin, no DownValues, not even
`Protected` — so it does nothing on its own; the graphics back end reads it out of a call's
option list.

`AspectRatio -> r` sets an explicit height-to-width ratio; `AspectRatio -> Automatic` takes
the true geometry from the coordinate values; `AspectRatio -> Full` stretches to fill the
region. `Plot` reports a stored default of `Automatic` and otherwise lays a standalone
curve out at `1/GoldenRatio`. The option changes the displayed shape only, never the
plotted data, so a call carrying it still evaluates to an ordinary `Graphics[...]` object.
