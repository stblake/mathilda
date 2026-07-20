# ListConvolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ListConvolve[ker, list]
    forms the convolution Sum_r ker[r] list[s-r] of ker with list.
ListConvolve[ker, list, k]
    aligns element k of ker with each element of list (cyclic).
ListConvolve[ker, list, {kL, kR}]
    sets the overhang: {-1,1} none (default), {1,1}/{-1,-1} maximal at one
    end, {1,-1} maximal at both.
ListConvolve[ker, list, klist, padding]
    pads list at each end with p, cyclic repetitions of {p1,p2,...}, the
    list itself (default), or {} for no padding.
ListConvolve[ker, list, klist, padding, g, h]
    uses g in place of Times and h in place of Plus.
ListConvolve[ker, list, klist, padding, g, h, lev]
    works at level lev. ker and list may be multidimensional. Large numeric
    input uses an FFT; exact and symbolic input is computed directly.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ListConvolve[{x, y}, {a, b, c, d, e, f}]
Out[1]= {b x + a y, c x + b y, d x + c y, e x + d y, f x + e y}

In[2]:= ListConvolve[{x, y, z}, {1, 2, 3, 4, 5, 6}, {1, -1}]
Out[2]= {x + 6 y + 5 z, 2 x + y + 6 z, 3 x + 2 y + z, 4 x + 3 y + 2 z, 5 x + 4 y + 3 z, 6 x + 5 y + 4 z, x + 6 y + 5 z, 2 x + y + 6 z}

In[3]:= ListConvolve[{{1, 1}, {1, 1}}, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 1]
Out[3]= {{20, 18, 22}, {14, 12, 16}, {26, 24, 28}}

In[4]:= ListConvolve[{x, y, z}, {1, 2, 3}, 1, {}, f, g]
Out[4]= {g[f[z], f[y], f[x, 1]], g[f[z], f[y, 1], f[x, 2]], g[f[z, 1], f[y, 2], f[x, 3]]}

In[5]:= ListCorrelate[{x, y}, {a, b, c, d, e, f}]
Out[5]= {a x + b y, b x + c y, c x + d y, d x + e y, e x + f y}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
