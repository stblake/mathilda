### Worked examples

```mathematica
In[1]:= Kurtosis[{1, 2, 3, 4, 5}]  (* the Pearson kurtosis of an evenly spaced integer sample *)
Out[1]= 17/10
```

```mathematica
In[1]:= Kurtosis[{1, 2, 4, 8}]  (* exact input gives an exact rational *)
Out[1]= 25141/13225
```

```mathematica
In[1]:= Kurtosis[{2., 4., 4., 4., 5., 5., 7., 9.}]  (* a machine-real sample gives a machine-real coefficient *)
Out[1]= 2.78125
```

### Notes

`Kurtosis[data]` is the standardized fourth central moment,
`CentralMoment[data, 4] / CentralMoment[data, 2]^2`. This is **Pearson** kurtosis,
not the excess form — subtract `3` to obtain the excess kurtosis (which is `0` for a
normal distribution). The exact value is a `Rational`, since both moments enter with
even powers.

For a matrix the coefficient is taken columnwise. A real `NDArray` / packed buffer
takes the kernel `ndred_kurtosis` and the head lowers inside `Compile[]`; an integer
buffer degrades to the exact `List` path.
