### Worked examples

```mathematica
In[1]:= FactorialPower[10, 3]
Out[1]= 720
```

```mathematica
In[1]:= FactorialPower[x, 4]
Out[1]= x (-3 + x) (-2 + x) (-1 + x)
```

```mathematica
In[1]:= FactorialPower[10, 4] == 10!/6!
Out[1]= True
```

### Notes

`FactorialPower[n, k]` is the falling factorial `n (n-1) (n-2) ... (n-k+1)`, the product of `k` descending linear factors. For symbolic `n` it expands to that product, e.g. `FactorialPower[x, 4] = x (x-1) (x-2) (x-3)`, the basis polynomial that makes the forward difference operator behave like differentiation (Δ x^(k) = k x^(k-1)). For non-negative integers it equals the ratio of factorials `n! / (n-k)!`, as the identity check confirms.
