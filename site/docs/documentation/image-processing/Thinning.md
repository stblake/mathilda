# Thinning

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Thinning[image] reduces the foreground to a one-pixel-wide skeleton by Zhang-Suen thinning, iterating until a pass deletes nothing. Thinning[image, n] stops after n iterations. The two subiterations are what preserve connectivity: deleting every individually-removable pixel in one pass severs a diagonal line, since two diagonal neighbours can each be removable while removing both disconnects the shape. A non-binary image is thresholded at 0.5 -- apply Binarize first for any other rule. The result is a "Bit" image, and it is always a subset of the input.`**

## Examples (26)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (6)

```mathematica
In[1]:= bar = Image[Table[If[3 <= i <= 7, 1, 0], {i, 1, 9}, {j, 1, 12}]];

In[2]:= Thinning[bar]
Out[2]= -Image-

In[3]:= {Round[Total[Flatten[ImageData[bar]]]], Round[Total[Flatten[ImageData[Thinning[bar]]]]]}
Out[3]= {60, 7}

In[4]:= ImageType[Thinning[bar]]
Out[4]= "Bit"

In[5]:= disk = Image[Table[N[Boole[(i - 12)^2 + (j - 12)^2 <= 81]], {i, 1, 24}, {j, 1, 24}]];

In[6]:= Thinning[disk]
Out[6]= -Image-
```

### Scope (8)

```mathematica
In[7]:= bar = Image[Table[If[3 <= i <= 7, 1, 0], {i, 1, 9}, {j, 1, 12}]];

In[8]:= diag = Image[Table[If[Abs[i - j] <= 1, 1, 0], {i, 1, 12}, {j, 1, 12}]];

In[9]:= Thinning[diag]
Out[9]= -Image-
```

An iteration limit stops early, so more foreground survives than at convergence

```mathematica
In[10]:= Round[Total[Flatten[ImageData[Thinning[bar, 1]]]]] >= Round[Total[Flatten[ImageData[Thinning[bar]]]]]
Out[10]= True
```

```mathematica
In[11]:= cross = Image[Table[If[Abs[i - 12] <= 2 || Abs[j - 12] <= 2, 1, 0], {i, 1, 24}, {j, 1, 24}]];

In[12]:= Thinning[cross]
Out[12]= -Image-

In[13]:= ring = Image[Table[N[Boole[36 <= (i - 12)^2 + (j - 12)^2 <= 100]], {i, 1, 24}, {j, 1, 24}]];

In[14]:= Thinning[ring]
Out[14]= -Image-
```

### Applications (4)

```mathematica
In[15]:= disk = Image[Table[N[Boole[(i - 12)^2 + (j - 12)^2 <= 81]], {i, 1, 24}, {j, 1, 24}]];
```

Skeleton then prune: the usual pairing, since a skeleton grows short spurs at boundary irregularities

```mathematica
In[16]:= Pruning[Thinning[disk], 2]
Out[16]= -Image-
```

How many branch pixels a shape's skeleton has

```mathematica
In[17]:= Round[Total[Flatten[ImageData[Thinning[disk]]]]]
Out[17]= 1
```

The skeleton of a binarised gradient, end to end

```mathematica
In[18]:= Thinning[Binarize[Image[Table[N[Boole[Abs[i - j] <= 3]], {i, 1, 20}, {j, 1, 20}], "Real"]]]
Out[18]= -Image-
```

### Properties & Relations (8)

```mathematica
In[19]:= bar = Image[Table[If[3 <= i <= 7, 1, 0], {i, 1, 9}, {j, 1, 12}]];

In[20]:= line = Image[Table[If[i == 5, 1, 0], {i, 1, 9}, {j, 1, 12}]];

In[21]:= diag = Image[Table[If[Abs[i - j] <= 1, 1, 0], {i, 1, 12}, {j, 1, 12}]];
```

Only ever deletes: the skeleton is a subset of what it came from

```mathematica
In[22]:= Max[Flatten[ImageData[Thinning[bar]] - ImageData[bar]]] <= 0
Out[22]= True
```

Settled means settled

```mathematica
In[23]:= ImageData[Thinning[Thinning[bar]]] === ImageData[Thinning[bar]]
Out[23]= True
```

Something already thin is untouched

```mathematica
In[24]:= ImageData[Thinning[line]] === ImageData[line]
Out[24]= True
```

THE DIAGONAL: still one connected component, which a single-pass thinning would break

```mathematica
In[25]:= Max[Flatten[MorphologicalComponents[Thinning[diag]]]]
Out[25]= 1
```

An isolated pixel has no neighbours, so no rule may remove it

```mathematica
In[26]:= Round[Total[Flatten[ImageData[Thinning[Image[Table[If[i == 5 && j == 6, 1, 0], {i, 1, 9}, {j, 1, 12}]]]]]]]
Out[26]= 1
```

## Algorithm

imagethin.c -- Thinning and Pruning: reducing a shape to its skeleton, and tidying it.

Both are ITERATIVE and both are defined on a BINARY image, which is what separates them from everything in imagefilter.c: a filter reads a window and writes one value, where these two delete pixels in passes and stop when a pass changes nothing. A neighbourhood kernel cannot express that — the result of one pass is the input to the next.

WHY ZHANG-SUEN. It is the standard two-subiteration thinning, and its two subiterations exist for a reason worth stating: deleting every deletable pixel in ONE pass severs a diagonal line, because two diagonal neighbours can each be individually removable while removing both disconnects the shape. Alternating the two conditions removes from opposite sides on alternating passes, which is what preserves connectivity. A single-pass "delete if removable" thinning looks correct on a thick blob and quietly breaks every diagonal stroke.

THRESHOLD. A non-binary image is thresholded at 0.5 rather than at "nonzero". Nonzero is the right rule for MorphologicalComponents, where the caller has usually binarised already, but for these two it would make almost every grey image entirely foreground and the skeleton would be a frame around the border. Callers wanting another rule should apply Binarize first, which is a decision they can see.

## Implementation notes

- `Protected`. Returns a `"Bit"` image, always a **subset** of the input.
- Zhang-Suen thinning: two subiterations per pass, deleting from opposite sides on alternate
  passes. That is what preserves **connectivity** — deleting every individually-removable
  pixel in a single pass severs a diagonal line, since two diagonal neighbours can each be
  removable while removing both disconnects the shape. A single-pass version looks correct on
  a thick blob and quietly breaks every diagonal stroke.
- A non-binary image is thresholded at `0.5`. Apply `Binarize` first for any other rule.
- Something already one pixel wide is left exactly alone, and an isolated pixel — no
  neighbours at all — is never removed.

**Attributes:** `Protected`.

## References

**See also:** [Binarize](../../image-processing/Binarize/)

- Source: [`src/imagethin.c`](https://github.com/stblake/mathilda/blob/main/src/imagethin.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
