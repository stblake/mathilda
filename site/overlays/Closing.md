### Worked examples

```mathematica
(* closing fills a single dark hole smaller than the element *)
In[1]:= ImageData[Closing[Image[{{1, 1, 1}, {1, 0, 1}, {1, 1, 1}}, "Bit"], 1]]
```

```mathematica
(* a "Bit" image closes to a "Bit" image *)
In[1]:= ImageType[Closing[Image[{{1, 1, 1}, {1, 0, 1}, {1, 1, 1}}, "Bit"], 1]]
```

```mathematica
(* idempotent: closing twice equals closing once *)
In[1]:= ImageData[Closing[Closing[Image[{{1, 1, 1}, {1, 0, 1}, {1, 1, 1}}, "Bit"], 1], 1]]
```

### Notes

`Closing[image, r]` dilates then erodes with the same element, filling dark
features smaller than it. Using the **same** structuring element for both passes
is what makes it idempotent — `Closing[Closing[f]] = Closing[f]` — and what makes
it a closing rather than two unrelated smoothings.

Like `Opening`, it is idempotent, and the two bracket the image:
`Erosion <= Opening <= image <= Closing <= Dilation` pointwise everywhere. Each
pass reuses `Dilation`/`Erosion`'s separable van Herk–Gil-Werman machinery, so
the cost does not grow with the radius. A `"Bit"` image stays `"Bit"`; other
types give `"Real"`.
