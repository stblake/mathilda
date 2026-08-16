# BoxMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BoxMatrix[r] gives a (2r+1) x (2r+1) matrix of 1s. It is NOT normalised, matching Mathematica, so ImageConvolve[image, BoxMatrix[1]] is nine times too bright; the normalised version is a mean filter. Kept faithful rather than helpfully rescaled, since a caller using BoxMatrix in arithmetic needs the ones.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (1)

```mathematica
In[1]:= BoxMatrix[1]
Out[1]= {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
