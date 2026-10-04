### Worked examples

```mathematica
(* a lone pixel grows to fill the 3x3 square neighbourhood *)
In[1]:= ImageData[Dilation[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}], 1]]
```

```mathematica
(* a structuring element uses its SUPPORT -- a cross grows a plus shape *)
In[1]:= ImageData[Dilation[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}, "Bit"], {{0, 1, 0}, {1, 1, 1}, {0, 1, 0}}]]
```

```mathematica
(* a "Bit" image dilates to a "Bit" image -- a max of 0/1 values is 0/1 *)
In[1]:= ImageType[Dilation[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}, "Bit"], 1]]
```

### Notes

`Dilation[image, r]` gives the maximum over a `(2r+1) × (2r+1)` square;
`Dilation[image, elem]` uses the support (nonzero positions) of the matrix
`elem`. This is **flat** morphology — the element's values do not enter the
maximum — which keeps `Dilation[img, BoxMatrix[1]]` and `Dilation[img, 1]` the
same operation.

A full rectangle is separable for the maximum, and each 1-D pass runs in three
comparisons per pixel by the van Herk–Gil-Werman algorithm, so the cost does not
grow with `r`. Padding replicates the border, so `Dilation >= image` holds at the
edges too. A `"Bit"` image gives a `"Bit"` image; other types give `"Real"`.
Together with `Erosion`, `Opening` and `Closing` it satisfies
`Erosion <= Opening <= image <= Closing <= Dilation` pointwise everywhere.
