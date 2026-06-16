### Worked examples

```mathematica
In[1]:= RepeatedTiming[Sum[i, {i, 100}]]
Out[1]= {4.53806e-05, 5050}
```

```mathematica
In[1]:= Last[RepeatedTiming[Factorial[100]]]
Out[1]= 93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000
```

```mathematica
In[1]:= Last[RepeatedTiming[Total[Range[10000]]]]
Out[1]= 50005000
```

### Notes

`RepeatedTiming` evaluates the expression many times and returns `{averageSeconds, result}`, giving a steadier estimate than `Timing` for fast operations. The timing value varies between runs; only the result element is reproducible.
