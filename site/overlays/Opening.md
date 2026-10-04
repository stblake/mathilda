### Worked examples

```mathematica
In[1]:= ImageData[Opening[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}], 1]]  (* a lone bright pixel is smaller than the element, so it is removed *)
```

```mathematica
In[1]:= blob = Image[{{0, 0, 0, 0}, {0, 1, 1, 0}, {0, 1, 1, 0}, {0, 0, 0, 0}}];  (* a 2x2 bright blob *)
```

```mathematica
In[1]:= ImageData[Opening[blob, 0]] == ImageData[Opening[Opening[blob, 0], 0]]  (* idempotent: opening twice is opening once *)
```

### Notes

`Opening[image, r]` erodes then dilates with the same element, removing bright features smaller
than it while leaving larger ones close to their original size.

It is **idempotent**: `Opening[Opening[f]] == Opening[f]`, the defining property of an opening
and the reason opening twice is not a sharpening loop — which is why the same element must be
used for both passes. Opening brackets the image from below in the morphology ordering
`Erosion ≤ Opening ≤ image ≤ Closing ≤ Dilation`. A `"Bit"` image stays `"Bit"`; other types
give `"Real"`.
