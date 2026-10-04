# ImageSize

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`ImageSize`**

is an option for Graphics and Plot that specifies the overall size of the image to display.

<details>
<summary>Notes</summary>

ImageSize -\> w sets the width to w pixels, with the height following from AspectRatio; ImageSize -\> {w, h} fixes both the width and height in pixels. The default is a width of 800.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A bare option symbol, not a function

```mathematica
In[1]:= Head[ImageSize]
Out[1]= Symbol
```

It only ever sits on the left of a rule

```mathematica
In[2]:= FullForm[ImageSize -> 400]
Out[2]= Rule[ImageSize, 400]
```

Accepted, still a Graphics

```mathematica
In[3]:= Head[Graphics[{Disk[]}, ImageSize -> {600, 400}]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `ImageSize` is an **option name** for `Graphics` and `Plot` giving the
overall size of the image to display. It is a bare option symbol, not a function: no
builtin, no DownValues, and not even the `Protected` attribute — its only registration is
the docstring set in `info_init` (`src/info.c`). Its meaning lives in the graphics back end
(`src/graphics/`), which reads it from a call's option list when sizing the window or
export.

**Representation.** `ImageSize` stays an inert `EXPR_SYMBOL` and appears only on the left
of a rule, so `FullForm[ImageSize -> 400]` is `Rule[ImageSize, 400]`. The value is either a
width `w` in pixels (the height then following from `AspectRatio`) or a `{w, h}` pair
fixing both. The default width is 800.

**Usage & limits.** Supplied as `Graphics[prims, ImageSize -> 400]` or
`Plot[f, {x, a, b}, ImageSize -> {600, 400}]`, it is consumed only at render/export time
and never changes the plotted data — a call carrying it still reduces to an ordinary
`Graphics[...]` object. It is not among the rules returned by `Options[Plot]`, since `Plot`
forwards it to the renderer rather than storing it as a plotting default.

**Attributes:** none registered.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`ImageSize` is an option for `Graphics` and `Plot` that sets the overall display size. It
is an inert option-name symbol — no builtin, no DownValues, not even `Protected` — so it
does nothing on its own; the graphics back end reads it from a call's option list when it
sizes the window or an export.

`ImageSize -> w` sets the width in pixels, with the height following from `AspectRatio`;
`ImageSize -> {w, h}` fixes both. The default width is 800. The option affects only the
rendered size, never the plotted data, so a call carrying it still evaluates to an ordinary
`Graphics[...]` object.
