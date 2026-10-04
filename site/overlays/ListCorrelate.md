### Worked examples

```mathematica
In[1]:= ListCorrelate[{1, 1}, {1, 2, 3, 4}]  (* no overhang: sliding pair sums *)
```

```mathematica
In[1]:= ListCorrelate[{1, 1, 1}/3, {1, 2, 3, 4, 5}]  (* a normalised kernel is a moving average *)
```

```mathematica
In[1]:= ListCorrelate[{1, 1}, {1, 2, 3, 4}, 1]  (* a single index argument makes the window cyclic *)
```

```mathematica
In[1]:= ListCorrelate[{x, y}, {a, b, c}]  (* symbolic data, aligned list[s+r] rather than list[s-r] *)
```

### Notes

`ListCorrelate[ker, list]` forms `Sum_r ker[r] list[s+r]`. It differs from `ListConvolve` only
in the sign of the kernel offset, so for one-dimensional data `ListCorrelate[ker, list]`
equals `ListConvolve[Reverse[ker], list]`. With no overhang (the default `{1, -1}`) the output
has length `Length[list] - Length[ker] + 1`; a normalised constant kernel therefore gives a
moving average.

The argument forms are identical to `ListConvolve` — a single integer `k` (cyclic alignment of
element `k`), a `{kL, kR}` overhang, a padding argument (constant, cyclic pad list, the list
itself by default, or `{}` for none), and the optional `g`/`h` operators and level. The
overhang settings are *negated* relative to `ListConvolve`: `{1, -1}` is none here and `{-1, 1}`
is maximal at both ends.

Data may be symbolic, exact, machine, or arbitrary-precision, and may be multidimensional.
Exact input stays exact; large numeric input with the default `Times`/`Plus` uses an FFT; a
packed kernel or list is read from its buffer and gives a packed result.
