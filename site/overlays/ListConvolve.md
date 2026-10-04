### Worked examples

```mathematica
In[1]:= ListConvolve[{1, 1}, {1, 2, 3, 4}]  (* no overhang: output is shorter than the list *)
```

```mathematica
In[1]:= ListConvolve[{1, -1}, {1, 2, 3, 4}]  (* a difference kernel gives successive differences *)
```

```mathematica
In[1]:= ListConvolve[{x, y}, {a, b, c}]  (* symbolic data is handled by the direct engine *)
```

```mathematica
In[1]:= ListConvolve[{1, 1, 1}, {1, 2, 3, 4, 5}, {1, -1}]  (* maximal overhang at both ends, cyclic wrap *)
```

### Notes

`ListConvolve[ker, list]` forms `Sum_r ker[r] list[s-r]` over the alignment window set by the
overhang. With no overhang (the default `{-1, 1}`) the output has length
`Length[list] - Length[ker] + 1`.

The third argument sets the overhang. A single integer `k` aligns element `k` of the kernel
with each list element cyclically; `{kL, kR}` places `ker[[kL]]` against the first output and
`ker[[kR]]` against the last. Common settings: `{-1, 1}` none (default), `{1, 1}`/`{-1, -1}`
maximal at one end, `{1, -1}` maximal at both. A fourth argument gives the padding — a constant,
a cyclic pad list, the list itself (the default, so out-of-range indices wrap), or `{}` for no
padding (the missing list factor is simply dropped). A fifth and sixth replace `Times`/`Plus`
with your own `g`/`h`, and a seventh sets the level.

Data may be symbolic, exact, machine, or arbitrary-precision, and kernel and list may be
multidimensional (the overhang broadcasts over axes, or `{{kL…},{kR…}}` sets each axis). Exact
input stays exact. Large numeric input with the default `Times`/`Plus` is computed by an FFT;
packed kernels and lists are read from their buffers and give a packed result. For
one-dimensional data `ListCorrelate[ker, list]` equals `ListConvolve[Reverse[ker], list]`.
