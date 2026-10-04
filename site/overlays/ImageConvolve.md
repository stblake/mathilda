### Worked examples

```mathematica
(* the identity kernel leaves a delta image unchanged *)
In[1]:= ImageData[ImageConvolve[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}]]
```

```mathematica
(* convolution reflects the kernel: a delta with {{1,2,3}} gives {1,2,3} ... *)
In[1]:= ImageData[ImageConvolve[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], {{1, 2, 3}}]]
```

```mathematica
(* ... where ImageCorrelate, which does not reflect, gives {3,2,1} *)
In[1]:= ImageData[ImageCorrelate[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], {{1, 2, 3}}]]
```

```mathematica
(* a constant image through a kernel summing to 1 is unchanged, border included *)
In[1]:= ImageData[ImageConvolve[Image[{{0.5, 0.5}, {0.5, 0.5}}], GaussianMatrix[1]]]
```

### Notes

`ImageConvolve[image, kernel]` convolves with a rank-2 numeric kernel. It is
**true convolution**: the kernel is reflected on both axes before summing, so it
differs from `ImageCorrelate` on an asymmetric kernel (the two agree exactly on a
symmetric one such as a Gaussian or a box).

Out-of-range reads clamp to the nearest edge pixel (`"Fixed"` padding), so a
constant image convolved with a kernel summing to 1 comes back unchanged
everywhere, including the edges — zero padding would darken them. The result is
always a `"Real"` image of the same dimensions, since a filtered byte is not
generally a byte, and each colour channel is convolved independently. A separable
kernel is factorised automatically into two 1-D passes; a volume takes a rank-3
convolution, and a rank-3 kernel handed to a plane is declined rather than
guessed at.
