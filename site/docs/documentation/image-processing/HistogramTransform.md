# HistogramTransform

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HistogramTransform[image] equalises the histogram, spreading the brightness distribution toward uniform over 256 bins by mapping each value through the cumulative distribution. The mapping is computed from the LUMINANCE and applied to every channel as a ratio, so hue survives; equalising each channel independently would shift colour, since it removes exactly the imbalance that makes an image warm or cool. A black pixel has no ratio to scale and takes the new luminance in every channel. Alpha passes through.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= dark = Image[Table[N[(i + j)/64], {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[2]:= HistogramTransform[dark]
Out[2]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABdElEQVR42u3dOW7DUAwAUX8t1i7d/7TJBcgiKfgd5E3JwnA0IDFADLu11r5eAcMwROPXOI7hfJqmcD7Pczh/v9/hfFmWcL6uazjfti2c7/sezo/jCOfneYbz67rC+X3f4fx5nh/N46eMMggggAAQ8H+Z1E5N7aggJwgEEAACPq2C1E5N7WSvbwOcIAJAAAHoVUFqp6Z2svdjA5wgAkAAAehVQWqnpnay928DnCACQAAB6FVBaqemdrK/1wY4QQSAAALQq4LUTk3tZM/HBjhBBIAAAtCrgtROTe1kz9MGOEEEgAAC0KuC1E5N7WS1aQOcIAJAAAHoVUFqp6Z2sv882gAniAAQQAB6VZDaqamd7FPoNsAJIgAEEIBeFaR2amon+0YCG+AEEQACCECvClI7NbWTfTulDXCCCAABBKBXBamdmtrJfqnEBjhBBIAAAtCrgtROUe0kv1prA5wgAkAAAfgrFaR2flc7rTUb4ASBAAJAwCfxDbAPIb36YSl0AAAAAElFTkSuQmCC)

Nothing in the original is above 0.5; afterwards the range is covered

```mathematica
In[3]:= {Round[Max[Flatten[ImageData[dark]]], 0.01], Round[Max[Flatten[ImageData[HistogramTransform[dark]]]], 0.01], Round[Min[Flatten[ImageData[HistogramTransform[dark]]]], 0.01]}
Out[3]= {0.5, 1.0, 0.0}
```

```mathematica
In[4]:= HistogramTransform[Image[Table[{N[i/32], N[j/32], 0.25}, {i, 1, 16}, {j, 1, 16}], "Real"]]
Out[4]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAADqElEQVR42u2dSY7bQBAEe1+4ShrY44ON+f8PDXiG9AcqDg0YEAVnHgsUtSRZCmQv9N6n0xnyoTuzXma7XqtZD9vdru+rWXdrto//cbPrP+2634Nd/7Df1/+yv5e7nfb7fkz2ed7t3+F0cB4nPVUyQAbIAEkG/L9K3iXbGW/TSAxQzzYVhNP+90/Jft80A6V4bx8f7Gso7osDHLGPh88ZNqCdbv8OHmmn6g5QC5JkgAyQZMClKCgC7aRQbFrwtmcl2v/yuS/gvIfzRLu+7/bn6Q0oDqhpt6kmThWoBq7cFOB7JXjFoTtALUiSATJAkgGXoqDkbAoqQEclTlC3qSlDdtQmO/MpfRrKfEqGTKkD1UAWlDNkSg2o5oCMC65p74LuALUgSQbIAEkGXIqCirfppQabIlqkun2eVhrQkZ35TNk+T2kNhvQgU0oR6KgMXYkJzhP7WOYTXdQdoBYkyQAZIMmAS1FQhcynUx1GyhpkPjPQ0dQgUwKq6TCPqDSYbwOZTyU6qjDPBzKfHCnbCTAd6Ut3gFqQJANkgCQDLkVB3dmUMkMWRFQzpwZ0ZFPHmuzz9GqvTctERzCPKEPmE+GHoPlIiTKfrwMmWoWha113gFqQDJBkgAyQnkVBE4yIzZAFLVDfgI7WYmc+Fahmhcyn0Up8WJM1ER1B5kNr2SplPlA/YETMwzxr3QFqQTJAkgEyQHoWBS3eposdsqAt2pnPBJnPAzKfubahEbEZqKbWArOUHcxfsq+53GhngFF6CdovSC1IkgEyQJIBL0FBG2Y+NkfcIfO5wX5BHajmFmHeEcyOppGvSivrSx4a+cq07xB8fsp2KAsSBakFSTJABkgy4GoUtMNufg/YO/oBGRFRzZrtzGeJCUayYEQM6IWynUJ1oJ1MmQ9tGARUc8DRX6IgtSBJBsgASQZcjILuMC/oAWvB3qB+hzViK1DNQivlgXZGs6CKdER7UPuhXRMPpCC7/qldE9WCJBkgAyQZcDEKegMK+uaJgvJQFrSFBPN80hgFQYbToF4x8wlDa8fc0K5Azn0CBf1RFqQWJMkAGSDJgFehoDBIQUA7axjMgsJYFkQjYpQF0fPIGHfOoRExynx+a78gtSBJBsgASQZcjILeYY3YdxwpG6OgDehlRgoKQ0/QcFQfpZ1R0RM0gk1Bi0bE1IIkGSADJBlwNQpylAVloCCbdm6wpmyFOu4uCHTk/hUF+aFoh7OgA8bEwqmnqaoFSTJABkgy4BX0FzFFT/9u36+kAAAAAElFTkSuQmCC)

### Properties & Relations (4)

```mathematica
In[5]:= dark = Image[Table[N[(i + j)/64], {i, 1, 16}, {j, 1, 16}], "Real"];
```

Monotone: equalisation is a cumulative distribution, so it cannot reorder pixels

```mathematica
In[6]:= Module[{d = ImageData[HistogramTransform[dark]]}, d[[1, 1]] <= d[[8, 8]] && d[[8, 8]] <= d[[16, 16]]]
Out[6]= True
```

```mathematica
In[7]:= ImageDimensions[HistogramTransform[dark]] === ImageDimensions[dark]
Out[7]= True
```

An already-spread image is near a fixed point, which is what "toward uniform" means

```mathematica
In[8]:= Module[{a = HistogramTransform[dark]}, Max[Abs[Flatten[ImageData[HistogramTransform[a]] - ImageData[a]]]] < 0.2]
Out[8]= True
```

## Algorithm

imagecolor.c -- ColorReplace, ColorQuantize and HistogramTransform.

Three heads that act on an image's COLOURS rather than its geometry, and they share the one thing that makes such operations awkward: a decision made per pixel needs a global view first. Replacing a colour needs a distance rule, quantising needs a palette derived from every pixel, and equalising needs the whole distribution. So each of these makes a pass to gather, then a pass to write — which is why none of them fits the filter machinery in imagefilter.c.

## Implementation notes

- `Protected`. Each value is mapped through the cumulative distribution over 256 bins, spreading
  the histogram toward uniform.
- The mapping is computed from the **luminance** and applied to every channel as a ratio, so hue
  survives. Equalising each channel independently would shift colour — it removes exactly the
  imbalance that makes an image warm or cool.
- A black pixel has no ratio to scale and takes the new luminance in every channel.
- **Monotone**: a cumulative distribution can never reorder two pixels.
- Alpha passes through.

**Attributes:** `Protected`.

## References

- Source: [`src/imagecolor.c`](https://github.com/stblake/mathilda/blob/main/src/imagecolor.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
