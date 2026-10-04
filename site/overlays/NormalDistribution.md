### Worked examples

```mathematica
In[1]:= NormalDistribution[2, 3]  (* a specified distribution prints its parameters in full *)
```

```mathematica
In[1]:= PDF[NormalDistribution[0, 1], 0]  (* the standard-normal density at 0 is 1/Sqrt[2 Pi] *)
```

```mathematica
In[1]:= PDF[NormalDistribution[0, 1], 1]  (* and away from the mean *)
```

```mathematica
In[1]:= PDF[NormalDistribution[], 0]  (* NormalDistribution[] is the standard normal *)
```

### Notes

`NormalDistribution[mu, sigma]` is a specification, not a value: it stays symbolic
and prints its parameters in full (unlike a fitted `LearnedDistribution`, which
prints elided). Pass it to `PDF[dist, x]` for the density `exp(-z^2/2)/(sigma
Sqrt[2 Pi])` with `z = (x - mu)/sigma`, or to `RandomVariate[dist, n]` for samples
(from the same stream as `RandomReal`, so `SeedRandom` makes them reproducible).
`NormalDistribution[]` is the standard normal. A non-positive `sigma`, or a
non-numeric evaluation point, leaves the consuming call unevaluated rather than
producing NaNs.
