---
status: Partial
references:
  - "Petkovšek, Wilf & Zeilberger, \"A=B\" (A K Peters, 1996)."
  - "Graham, Knuth & Patashnik, \"Concrete Mathematics\", 2nd ed. (Addison-Wesley, 1994), ch. 2 & 6."
---
### Worked examples

```mathematica
In[1]:= Sum[k, {k, 1, 10}]
Out[1]= 55
```

```mathematica
In[1]:= Sum[k^2, {k, 1, n}]
Out[1]= 1/6 n (1 + n) (1 + 2 n)
```

Faulhaber's formula extends to high powers, here the fifth power:

```mathematica
In[1]:= Sum[k^5, {k, 1, n}]
Out[1]= 1/12 n^2 (1 + n)^2 (-1 + 2 n + 2 n^2)
```

A finite geometric sum is returned in closed form in the parameter `r`:

```mathematica
In[1]:= Sum[r^k, {k, 0, n}]
Out[1]= -1/(-1 + r) + r^(1 + n)/(-1 + r)
```

Gosper's algorithm handles the arithmetic–geometric summand `k x^k`:

```mathematica
In[1]:= Sum[k x^k, {k, 1, n}]
Out[1]= x/(1 - 2 x + x^2) + (x^(1 + n) (-1 - n - x + (1 + n) x))/(1 - 2 x + x^2)
```

Convergent infinite geometric and exponential series close in symbolic form:

```mathematica
In[1]:= Sum[1/2^k, {k, 0, Infinity}]
Out[1]= 2

In[2]:= Sum[x^k/k!, {k, 0, Infinity}]
Out[2]= E^x
```

### Notes

`Sum` evaluates numeric ranges directly and closes symbolic finite ranges in closed form through the polynomial, geometric, and Gosper (`Method`) routines — so `Sum[k^2, {k, 1, n}]` returns Faulhaber's polynomial and `Sum[k x^k, {k, 1, n}]` is summed by the Gosper backend over a symbolic upper bound. Some infinite sums are recognised: geometric series such as `Sum[1/2^k, {k, 0, Infinity}]` give `2`, and the exponential generating function `Sum[x^k/k!, {k, 0, Infinity}]` returns `E^x`. Zeta-type series such as `Sum[1/k^2, {k, 1, Infinity}]` are **not** evaluated and stay symbolic. `Sum` is `HoldAll`, so the iterator variable is not evaluated before the range is set up.
