### Worked examples

```mathematica
In[1]:= MultiplicativeOrder[2, 7]
Out[1]= 3
```

```mathematica
In[1]:= MultiplicativeOrder[10, 7]
Out[1]= 6
```

```mathematica
In[1]:= MultiplicativeOrder[7, 1000000007]
Out[1]= 500000003
```

```mathematica
In[1]:= MultiplicativeOrder[3, 998244353]
Out[1]= 998244352
```

```mathematica
In[1]:= MultiplicativeOrder[2, 11, {1, 10}]
Out[1]= 5
```

### Notes

`MultiplicativeOrder[k, n]` is the smallest `m > 0` with `k^m ≡ 1 (mod n)`. The
order `6` for `10` modulo `7` reflects that `1/7 = 0.142857...` has a repeating
block of length `6`. The two large-modulus cases use prime moduli: `3` is a
primitive root of the NTT prime `998244353`, so its order equals `n - 1`. The
three-argument form `MultiplicativeOrder[k, n, {r1, ...}]` finds the least `m`
with `k^m` congruent to one of the listed residues. All arithmetic is exact via
GMP. The result is unevaluated when `gcd(k, n) ≠ 1`.
