# Pruning

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Pruning[image] removes one pixel from every free end of the foreground; Pruning[image, n] repeats that n times, which shortens each branch by up to n and deletes any branch shorter than that. Used after Thinning to remove the short spurs a skeleton grows at boundary irregularities. An end point has exactly one foreground neighbour, so an ISOLATED pixel is not one and survives: pruning shortens branches rather than erasing specks. Pruning[image, 0] is the image unchanged. The result is a "Bit" image.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= line = Image[Table[If[i == 5, 1, 0], {i, 1, 9}, {j, 1, 12}]];

In[2]:= {Round[Total[Flatten[ImageData[line]]]], Round[Total[Flatten[ImageData[Pruning[line]]]]]}
Out[2]= {12, 10}

In[3]:= Round[Total[Flatten[ImageData[Pruning[line, 3]]]]]
Out[3]= 6

In[4]:= ImageType[Pruning[line]]
Out[4]= "Bit"
```

### Properties & Relations (6)

```mathematica
In[5]:= line = Image[Table[If[i == 5, 1, 0], {i, 1, 9}, {j, 1, 12}]];

In[6]:= dot = Image[Table[If[i == 5 && j == 6, 1, 0], {i, 1, 9}, {j, 1, 12}]];
```

Zero passes is the identity -- the boundary a loop written with <= gets wrong

```mathematica
In[7]:= ImageData[Pruning[line, 0]] === ImageData[line]
Out[7]= True
```

An isolated pixel has no neighbours, so it is not an end point and survives any n

```mathematica
In[8]:= Round[Total[Flatten[ImageData[Pruning[dot, 5]]]]]
Out[8]= 1
```

Only ever deletes

```mathematica
In[9]:= Max[Flatten[ImageData[Pruning[line, 2]] - ImageData[line]]] <= 0
Out[9]= True
```

Each pass takes one pixel from each of the two ends

```mathematica
In[10]:= Table[Round[Total[Flatten[ImageData[Pruning[line, k]]]]], {k, 0, 4}]
Out[10]= {12, 10, 8, 6, 4}
```

## Algorithm

imagethin.c -- Thinning and Pruning: reducing a shape to its skeleton, and tidying it.

Both are ITERATIVE and both are defined on a BINARY image, which is what separates them from everything in imagefilter.c: a filter reads a window and writes one value, where these two delete pixels in passes and stop when a pass changes nothing. A neighbourhood kernel cannot express that — the result of one pass is the input to the next.

WHY ZHANG-SUEN. It is the standard two-subiteration thinning, and its two subiterations exist for a reason worth stating: deleting every deletable pixel in ONE pass severs a diagonal line, because two diagonal neighbours can each be individually removable while removing both disconnects the shape. Alternating the two conditions removes from opposite sides on alternating passes, which is what preserves connectivity. A single-pass "delete if removable" thinning looks correct on a thick blob and quietly breaks every diagonal stroke.

THRESHOLD. A non-binary image is thresholded at 0.5 rather than at "nonzero". Nonzero is the right rule for MorphologicalComponents, where the caller has usually binarised already, but for these two it would make almost every grey image entirely foreground and the skeleton would be a frame around the border. Callers wanting another rule should apply Binarize first, which is a decision they can see.

## Implementation notes

- `Protected`. Returns a `"Bit"` image, always a subset of the input.
- An end point has **exactly one** foreground neighbour, so an **isolated pixel is not one**
  and survives: pruning shortens branches rather than erasing specks. A rule that removed
  isolated pixels would quietly delete every one-pixel component.
- `n` passes shorten each branch by up to `n`, and remove any branch shorter than that
  entirely — which is what makes it the companion to `Thinning`, whose skeletons grow short
  spurs at boundary irregularities.
- `Pruning[image, 0]` is the image unchanged.

**Attributes:** `Protected`.

## References

**See also:** [Thinning](../../image-processing/Thinning/)

- Source: [`src/imagethin.c`](https://github.com/stblake/mathilda/blob/main/src/imagethin.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
