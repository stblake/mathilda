### Worked examples

```mathematica
(* equalisation preserves the image dimensions *)
In[1]:= ImageDimensions[HistogramTransform[Image[{{0.1, 0.5}, {0.5, 0.9}}]]]
```

```mathematica
(* the result is a "Real" image *)
In[1]:= ImageType[HistogramTransform[Image[{{0.1, 0.5}, {0.5, 0.9}}]]]
```

```mathematica
(* a one-channel grey image stays one channel *)
In[1]:= ImageChannels[HistogramTransform[Image[{{0.2, 0.4, 0.6, 0.8}}]]]
```

### Notes

`HistogramTransform[image]` equalises the histogram, spreading the brightness
distribution toward uniform over 256 bins by mapping each value through the
cumulative distribution.

The mapping is computed from the **luminance** and applied to every channel as a
ratio, so hue survives — equalising each channel independently would shift
colour, since it removes exactly the imbalance that makes an image warm or cool.
A black pixel has no ratio to scale and takes the new luminance in every channel.
An alpha channel passes through.
