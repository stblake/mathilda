### Worked examples

```mathematica
In[1]:= PrimePi[100]
Out[1]= 25
```

The count scales to large bounds; there are 78498 primes below one million:

```mathematica
In[1]:= PrimePi[10^6]
Out[1]= 78498
```

### Notes

`PrimePi[x]` gives the prime-counting function `π(x)`, the number of primes
less than or equal to `x`. For `x = 100` the answer is `25`; for `x = 10^6` it
is `78498`, consistent with the asymptotic `π(x) ~ x / Log[x]`.
