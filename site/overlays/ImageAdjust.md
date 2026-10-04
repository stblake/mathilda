### Worked examples

```mathematica
(* full-range stretch: the darkest pixel goes to 0, the brightest to 1 *)
In[1]:= ImageData[ImageAdjust[Image[{{0.2, 0.6}}]]]
```

```mathematica
(* a constant image has no range to stretch and is returned unchanged *)
In[1]:= ImageData[ImageAdjust[Image[{{0.5, 0.5}}]]]
```

```mathematica
(* parametric {contrast, brightness}: brightness is a plain offset about mid-grey *)
In[1]:= ImageData[ImageAdjust[Image[{{0.5}}], {0, 0.2}]]
```

### Notes

`ImageAdjust[image]` stretches to the full range: the darkest pixel becomes
exactly `0` and the brightest exactly `1`. It is idempotent — a second stretch is
the identity — and a constant image has no range to stretch and comes back
unchanged, since dividing by zero is not the answer.

`ImageAdjust[image, {c, b}]` and `[image, {c, b, g}]` apply contrast `c`,
brightness `b` and gamma `g` by the stated curve `v' = (v − 1/2)(1 + c) + 1/2 +
b`, clipped to `[0, 1]`, then raised to the power `1/g`. Contrast pivots about
mid-grey so it does not also shift brightness; clipping precedes gamma because a
negative base has no real power. This curve is Mathilda's documented choice, not
a claim of bit-compatibility with Mathematica. It accepts volumes as well as
planes.
