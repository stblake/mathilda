### Worked examples

```mathematica
In[1]:= Skewness[{1, 2, 3, 10}]  (* a right-skewed integer sample; the exact value is a radical *)
Out[1]= 18/25 Sqrt[2]
```

```mathematica
In[1]:= Skewness[{1, 2, 3, 4, 5}]  (* a symmetric sample has zero skewness *)
Out[1]= 0
```

```mathematica
In[1]:= Skewness[{2., 4., 4., 4., 5., 5., 7., 9.}]  (* a machine-real sample gives a machine-real coefficient *)
Out[1]= 0.65625
```

### Notes

`Skewness[data]` is the standardized third central moment,
`CentralMoment[data, 3] / CentralMoment[data, 2]^(3/2)` — a dimensionless measure
of asymmetry. A positive value signals a longer right tail, a negative value a
longer left tail, and a symmetric distribution gives exactly `0`. Exact input
yields exact output, which is a radical in general.

For a matrix the coefficient is taken columnwise (the `CentralMoment` ratio
threads). A real `NDArray` / packed buffer takes the kernel `ndred_skewness` and
the head lowers inside `Compile[]`; an integer sample's skewness is a radical, so
an integer buffer degrades to the exact `List` path.
