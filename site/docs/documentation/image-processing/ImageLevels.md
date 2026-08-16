# ImageLevels

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageLevels[image] gives {{level, count}, ...}: the histogram as DATA, not a plot -- use Histogram over the result for a picture. ImageLevels[image, n] uses n bins. Levels are on the same unit scale as ImageData, so a level can be compared against a pixel value without rescaling. A "Bit" image uses its 2 natural levels and "Byte" its 256, because those ARE the distinct values; a "Real" image has no natural set and is binned into 256 over [0, 1]. The counts sum to the pixel count exactly, every pixel landing in one bin. Accepts volumes as well as planes.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (2)

```mathematica
In[1]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[2]:= ImageLevels[bit]
Out[2]= {{0.0, 32}, {1.0, 32}}
```

### Scope (1)

```mathematica
In[3]:= ImageLevels[Image[{{0, 1}, {1, 0}}]]
Out[3]= {{0.0, 2}, {1.0, 2}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Histogram](../../graphics/Histogram/), [ImageData](../../image-processing/ImageData/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
