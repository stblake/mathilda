### Worked examples

```mathematica
In[1]:= NextPrime[100]
Out[1]= 101

In[2]:= NextPrime[1000]
Out[2]= 1009

In[3]:= NextPrime[10, 3]
Out[3]= 17
```

### Notes

`NextPrime[x]` gives the smallest prime greater than `x`. The two-argument form `NextPrime[x, k]` gives the k-th prime after `x`, so `NextPrime[10, 3]` steps past 11 and 13 to reach 17.
