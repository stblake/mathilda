### Worked examples

```mathematica
In[1]:= PrimitiveRootList[7]
Out[1]= {3, 5}
```

It works for twice-an-odd-prime-power moduli too:

```mathematica
In[1]:= PrimitiveRootList[18]
Out[1]= {5, 11}
```

When the group is non-cyclic (15 is neither `4`, an odd prime power, nor twice one), the list is empty:

```mathematica
In[1]:= PrimitiveRootList[15]
Out[1]= {}
```

The number of primitive roots of a prime `p` is `EulerPhi[EulerPhi[p]] = EulerPhi[p - 1]`; for `p = 101` this gives 40, matching the list length:

```mathematica
In[1]:= Length[PrimitiveRootList[101]]
Out[1]= 40

In[2]:= EulerPhi[EulerPhi[101]]
Out[2]= 40
```

### Notes

`PrimitiveRootList[n]` returns every primitive root of `n` in canonical
residues `{1, ..., n-1}`, sorted. The list is non-empty only when `n` is `2`,
`4`, an odd prime power `p^k`, or `2 p^k`; otherwise it is `{}`. When primitive
roots exist there are exactly `EulerPhi[EulerPhi[n]]` of them, since the cyclic
group of order `EulerPhi[n]` has that many generators.
