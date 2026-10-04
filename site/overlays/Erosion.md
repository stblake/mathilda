### Worked examples

```mathematica
In[1]:= ImageData[Erosion[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}], 1]]  (* a lone bright pixel is eroded away *)
```

```mathematica
In[1]:= ImageType[Erosion[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}], 1]]  (* a Bit image stays Bit *)
```

```mathematica
In[1]:= ImageData[Erosion[Image[{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}}], 1]]  (* a solid field survives: replicate padding *)
```

### Notes

`Erosion[image, r]` gives the minimum over a `(2r+1) × (2r+1)` square; `Erosion[image, elem]`
uses the nonzero support of the matrix `elem` as the neighbourhood (flat morphology — the
element's values do not enter the minimum). Padding replicates the border, the same rule the
convolutions use, which is what makes `Erosion <= image` hold at the edges too.

Erosion is dual to Dilation: for a symmetric element `Erosion[f, k] = 1 - Dilation[1 - f, k]`
exactly. A `"Bit"` image gives a `"Bit"` image; other types give `"Real"`. A full rectangle is
separable for the minimum, and a van Herk–Gil–Werman pass makes the cost independent of the
radius.
