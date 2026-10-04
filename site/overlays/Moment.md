### Worked examples

```mathematica
In[1]:= Moment[{2, 4, 4, 4, 5, 5, 7, 9}, 2]  (* the raw second moment of an integer sample is exact *)
Out[1]= 29
```

```mathematica
In[1]:= Moment[{1., 2., 3., 4.}, 3]  (* a machine-real sample gives a machine-real moment *)
Out[1]= 25.0
```

```mathematica
In[1]:= Moment[{1, 2, 3, 4}, 1]  (* the first raw moment is just the mean *)
Out[1]= 5/2
```

```mathematica
In[1]:= Moment[{1, 2, 3, 4}, 0]  (* and the zeroth raw moment is 1 for any data *)
Out[1]= 1
```

```mathematica
In[1]:= Moment[{{1, 2}, {3, 4}, {5, 6}}, 2]  (* on a matrix the moment is taken columnwise over the first axis *)
Out[1]= {35/3, 56/3}
```

```mathematica
In[1]:= Moment[{a, b, c}, 2]  (* symbolic data stays symbolic *)
Out[1]= 1/3 (a^2 + b^2 + c^2)
```

```mathematica
In[1]:= Moment[{{1, 2}, {3, 4}, {5, 6}}, {1, 2}]  (* a list order gives the multivariate mixed raw moment *)
Out[1]= 232/3
```

### Notes

The raw (power) moment `Moment[data, r]` is `CentralMoment[data, r]` without the
mean subtraction: `μ_r = (1/n) Σ xᵢ^r`. Because there is no mean to subtract,
`Mean[data^r]` threads correctly at every rank, so the matrix and higher-rank
columnwise forms need no special handling — the `Listable` `Power` does the work.

A real `NDArray` / packed buffer takes the kernel `ndred_moment`, and the head
lowers inside `Compile[]` for a real vector with integer order. An integer buffer
falls back to the exact `List` path, because an integer sample's raw moment is a
`Rational` no machine slot can hold.

`Moment` carries the `NHoldAll` attribute (as in Mathematica), so `N` does not
thread into a symbolic order or symbolic data.
