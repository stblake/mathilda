### Worked examples

```mathematica
In[1]:= ImageData[MedianFilter[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], 1]]  (* a lone bright pixel vanishes exactly *)
```

```mathematica
In[1]:= MedianFilter[Image[{{0.2, 0.9, 0.3}, {0.8, 1., 0.1}, {0.4, 0.7, 0.6}}], 1]  (* an outlier-rejecting smooth *)
```

```mathematica
In[1]:= ImageDimensions[MedianFilter[Image[{{0., 0.}, {0., 0.}}], 1]]  (* same size as the input *)
```

### Notes

`MedianFilter[image, r]` replaces each pixel with the median over a `(2r+1) × (2r+1)`
neighbourhood. Unlike a Gaussian it removes an isolated outlier **exactly** rather than
attenuating and smearing it, which is what makes it the filter for salt-and-pepper noise.

It is the one filter here that is **not** separable: a sum, a max and a min all decompose
because they ignore grouping, but a median depends on a value's rank within the whole window —
the median of the row-medians of `{{1,2,9},{3,4,5},{6,7,8}}` is 4 where the true median is 5.
For an even window the lower middle is taken rather than the average of the two, so the output
is always one of the inputs.
