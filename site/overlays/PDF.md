### Worked examples

```mathematica
In[1]:= PDF[NormalDistribution[], 0]  (* the standard normal at its peak *)
```

```mathematica
In[1]:= PDF[NormalDistribution[], {-1., 0., 1.}]  (* a scalar distribution threads over a list of points *)
```

```mathematica
In[1]:= PDF[UniformDistribution[{0., 2.}], {-1., 1., 3.}]  (* zero outside the support, flat inside *)
```

```mathematica
In[1]:= m = LearnDistribution[{{1., 2.}, {2., 3.}, {3., 5.}, {4., 4.}, {5., 7.}, {6., 8.}}]  (* a fitted multinormal *)
```

```mathematica
In[1]:= PDF[m, {{3.5, 4.8}, {50., 50.}}]  (* a matrix threads to one density per row *)
```

### Notes

`PDF[NormalDistribution[], 0]` is `1/Sqrt[2 Pi]` numerically, which is how the
closed-form density is verified. A uniform density is `1/(b − a)` on its closed support
and exactly zero outside.

The reading of a list argument depends on the distribution, and it has to. For a
**scalar** distribution a list of `x` threads to a list of densities — the shape a caller
plotting a density wants. For a **point** distribution (a fitted multinormal, mixture or
kernel estimate) the argument is itself a coordinate vector, so a list is *one*
observation, and only a *matrix* threads to one density per row. The second row above,
far out in the tail, comes back as a vanishingly small but non-zero density because the
point kernels work in log space.
