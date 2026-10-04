### Worked examples

```mathematica
(* a 3-pixel horizontal branch: both free ends are pruned, leaving the middle pixel *)
In[1]:= ImageData[Pruning[Image[{{0, 0, 0, 0}, {1, 1, 1, 0}, {0, 0, 0, 0}}], 1]]
```

```mathematica
(* Pruning[image, 0] returns the image unchanged *)
In[1]:= ImageData[Pruning[Image[{{0, 0, 0, 0}, {1, 1, 1, 0}, {0, 0, 0, 0}}], 0]]
```

```mathematica
(* the result is a "Bit" image *)
In[1]:= ImageType[Pruning[Image[{{0, 0, 0}, {1, 1, 1}, {0, 0, 0}}]]]
```

### Notes

`Pruning[image]` removes one pixel from every free end of the foreground;
`Pruning[image, n]` repeats that `n` times, shortening each branch by up to `n`
and deleting any branch shorter than that. `Pruning[image, 0]` is the image
unchanged.

It is used after `Thinning` to remove the short spurs a skeleton grows at
boundary irregularities. An end point has **exactly one** foreground neighbour,
so an isolated pixel (with none) is not one and survives: pruning shortens
branches rather than erasing specks. As in `Thinning`, each pass marks then
deletes together, so one removal does not change a neighbour's verdict mid-pass.
The result is a `"Bit"` image.
