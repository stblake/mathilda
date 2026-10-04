### Worked examples

```mathematica
In[1]:= CentralMoment[{2, 4, 4, 4, 5, 5, 7, 9}, 2]  (* the second central moment — Variance without the bias correction *)
Out[1]= 4
```

```mathematica
In[1]:= Variance[{2, 4, 4, 4, 5, 5, 7, 9}]  (* compare: Variance divides by n - 1, not n *)
Out[1]= 32/7
```

```mathematica
In[1]:= CentralMoment[{1., 2., 3., 4., 5.}, 2]  (* a machine-real sample gives a machine-real moment *)
Out[1]= 1.25
```

```mathematica
In[1]:= CentralMoment[{2, 4, 4, 4, 5, 5, 7, 9}, 1]  (* the first central moment is always zero *)
Out[1]= 0
```

```mathematica
In[1]:= CentralMoment[{1, 2, 3, 4}, 3]  (* an odd central moment vanishes on symmetric data *)
Out[1]= 0
```

```mathematica
In[1]:= CentralMoment[{{1, 2}, {3, 4}, {5, 6}}, 2]  (* on a matrix the moment is columnwise over the first axis *)
Out[1]= {8/3, 8/3}
```

```mathematica
In[1]:= CentralMoment[{{0, 1}, {2, 5}, {4, 3}}, {1, 1}]  (* a list order gives the multivariate mixed central moment *)
Out[1]= 4/3
```

### Notes

`CentralMoment[data, r] = (1/n) Σ (xᵢ − Mean[data])^r`. It is `Variance` with the
bias correction removed: it divides by `n` rather than `n − 1`, raises to the power
`r` rather than squaring, needs only `n ≥ 1`, and applies no `Conjugate`. So
`CentralMoment[data, 2]` and `Variance[data]` differ by exactly the factor
`(n − 1)/n`.

For a matrix or array the centering is done per first-axis slice, because the
obvious `data − Mean[data]` would thread row-wise rather than column-wise.

A real `NDArray` / packed buffer takes the kernel `ndred_central_moment` and the
head lowers inside `Compile[]`; an integer buffer degrades to the exact `List`
path, since the exact answer is a `Rational`.
