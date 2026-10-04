### Worked examples

```mathematica
In[1]:= k = SmoothKernelDistribution[Table[1. i/10., {i, 0, 40}]]  (* a KDE over 41 points on 0..4 *)
```

```mathematica
In[1]:= PDF[k, {2.0}]  (* evaluate the estimated density *)
```

```mathematica
In[1]:= k2 = SmoothKernelDistribution[Table[1. i/10., {i, 0, 40}], 0.5]  (* set the bandwidth explicitly *)
```

```mathematica
In[1]:= PDF[k2, {2.0}]  (* a different bandwidth gives a different estimate *)
```

### Notes

The sample *is* the model: nothing is fitted except the bandwidth, so a KDE costs
nothing to build and everything to evaluate — the same trade a nearest-neighbour
predictor makes. The kernel is a product Gaussian with a per-dimension bandwidth.

The default bandwidth is the multivariate normal-reference rule, which in one dimension
is exactly Silverman's `1.06 σ n^(−1/5)`. Being a normal-reference rule it oversmooths
strongly multimodal data — a known property of the rule rather than a defect — which is
why the explicit-bandwidth form shown above exists. A constant column has no scale, so
the rule gives a zero bandwidth and the call returns unevaluated rather than dividing by
zero.
