### Worked examples

```mathematica
In[1]:= Fibonacci[10]
Out[1]= 55
```

```mathematica
In[1]:= Fibonacci[100]
Out[1]= 354224848179261915075
```

```mathematica
In[1]:= Fibonacci[10, x]
Out[1]= 5 x + 20 x^3 + 21 x^5 + 8 x^7 + x^9
```

```mathematica
In[1]:= Fibonacci[200]/Fibonacci[199] // N
Out[1]= 1.61803
```

```mathematica
In[1]:= Sum[Fibonacci[k], {k, 1, 10}] == Fibonacci[12] - 1
Out[1]= True
```

### Notes

`Fibonacci[n]` uses GMP fast-doubling, so even `Fibonacci[100]` (21 digits) is returned exactly and instantly. The two-argument form gives the Fibonacci polynomial `F_n(x)` from the recurrence `F_k = x F_{k-1} + F_{k-2}`; `F_10(x)` factors the integer Fibonacci numbers (`F_10(1) = 55`). The ratio of consecutive terms converges to the golden ratio `φ = 1.61803...`, and the telescoping identity `∑_{k=1}^{n} F_k = F_{n+2} - 1` holds exactly.
