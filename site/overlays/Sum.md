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

```mathematica
In[1]:= Sum[k^3, {k, 1, n}]
Out[1]= 1/4 n^2 (1 + n)^2
```

```mathematica
In[1]:= Sum[2^k, {k, 0, 5}]
Out[1]= 63
```

### Notes

`Sum` evaluates numeric ranges directly and closes symbolic finite ranges in closed form through the polynomial, geometric, and Gosper (`Method`) routines — so `Sum[k^2, {k, 1, n}]` returns Faulhaber's polynomial. The Gosper backend handles hypergeometric summands over a symbolic upper bound. Infinite sums are **not** evaluated: `Sum[1/k^2, {k, 1, Infinity}]` and `Sum[1/2^k, {k, 0, Infinity}]` both return unevaluated, so convergent series like ζ(2) and geometric series stay symbolic. `Sum` is `HoldAll`, so the iterator variable is not evaluated before the range is set up.
