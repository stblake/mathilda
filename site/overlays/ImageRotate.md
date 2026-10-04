### Worked examples

```mathematica
In[1]:= ImageDimensions[ImageRotate[Image[{{0., 1., 0.}}]]]  (* a quarter turn swaps the dimensions: {3,1} -> {1,3} *)
```

```mathematica
In[1]:= ImageData[ImageRotate[Image[{{0., 1.}, {2., 3.}}]]]  (* a quarter turn counterclockwise, exact *)
```

```mathematica
In[1]:= ImageDimensions[ImageRotate[Image[{{0., 1., 0.}}], 90 Degree]]  (* an angle in degrees is accepted and numericalised *)
```

```mathematica
In[1]:= ImageData[ImageRotate[Image[{{0., 1.}, {2., 3.}}], Pi]]  (* a half turn is still an exact permutation *)
```

### Notes

`ImageRotate[image]` turns a quarter counterclockwise; `ImageRotate[image, angle]` turns by
`angle` radians (use `n Degree` for degrees, a negative angle turns clockwise); a side argument
turns the top of the image to face `Left`/`Right`/`Top`/`Bottom`, and `side1 -> side2` turns
one onto the other.

A multiple of a right angle takes an **exact** index-permutation path — nothing is
interpolated, four quarter turns are the identity, and an odd number of quarters swaps the
dimensions. Any other angle interpolates bilinearly, and area rotated in from outside reads as
0 rather than the replicated edge, since it was never photographed. Unlike Mathematica, a free
angle keeps the input's dimensions rather than enlarging to enclose the rotated image.
