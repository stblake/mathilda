### Worked examples

```mathematica
In[1]:= NextPrime[100]
Out[1]= 101

In[2]:= NextPrime[10, 3]
Out[2]= 17
```

The prime just above `10^20` — arbitrary-precision, so no overflow:

```mathematica
In[1]:= NextPrime[10^20]
Out[1]= 100000000000000000039
```

Stepping five primes past `2^31`:

```mathematica
In[1]:= NextPrime[2^31, 5]
Out[1]= 2147483777
```

### Notes

`NextPrime[x]` gives the smallest prime greater than `x`. The two-argument form
`NextPrime[x, k]` gives the k-th prime after `x`, so `NextPrime[10, 3]` steps past
11 and 13 to reach 17. Arguments may be arbitrarily large integers.
