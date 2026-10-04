### Worked examples

```mathematica
In[1]:= Head[UniformDistribution[{0, 1}]]  (* an inert distribution object *)
```

```mathematica
In[1]:= PDF[UniformDistribution[{0, 1}], 1/2]  (* density 1/(hi - lo) inside the support *)
```

```mathematica
In[1]:= PDF[UniformDistribution[{2, 6}], 10]  (* zero strictly outside [lo, hi] *)
```

### Notes

`UniformDistribution[{lo, hi}]` represents the continuous uniform distribution on
`[lo, hi]`, and `UniformDistribution[]` the standard uniform on `{0, 1}`. It is an inert,
`Protected` object — it evaluates to itself, `Head` is `UniformDistribution` — whose meaning
is supplied by the consumers `PDF` and `RandomVariate`.

`PDF[UniformDistribution[{lo, hi}], x]` is `1/(hi - lo)` for `lo <= x <= hi` (endpoints
included) and `0` outside; these numeric consumers work at machine precision, so the density
comes back as a machine real. Symbolic moment functions such as `Mean`, `Variance` and `CDF`
are not yet wired to this object and are left unevaluated.
