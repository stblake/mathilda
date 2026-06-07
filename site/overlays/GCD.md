---
status: Stable
references:
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", on the Euclidean algorithm."
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\", on GCD computation over the integers and rationals."
---
### Worked examples

```mathematica
In[1]:= GCD[12, 18, 30]
Out[1]= 6
```

```mathematica
In[1]:= GCD[2^20, 2^15]
Out[1]= 32768
```

```mathematica
In[1]:= GCD[1/2, 1/3]
Out[1]= 1/6
```

```mathematica
In[1]:= GCD[0, 5]
Out[1]= 5
```

### Notes

GCD folds the Euclidean algorithm across all arguments, so three-or-more-argument
calls such as `GCD[12, 18, 30]` reduce pairwise to `6`. It extends to rationals
via `gcd(a/b, c/d) = gcd(a,c)/lcm(b,d)`, giving `GCD[1/2, 1/3] = 1/6`. The
convention `GCD[0, n] = n` holds, since zero is divisible by every integer. Large
powers of two are handled exactly through GMP, with `GCD[2^20, 2^15] = 2^15 =
32768`.
