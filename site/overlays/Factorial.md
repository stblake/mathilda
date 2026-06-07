---
status: Stable
references:
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", on factorial computation."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on the Gamma function and special values."
---
### Worked examples

```mathematica
In[1]:= 20!
Out[1]= 2432902008176640000
```

```mathematica
In[1]:= Factorial[30]
Out[1]= 265252859812191058636308480000000
```

```mathematica
In[1]:= 0!
Out[1]= 1
```

```mathematica
In[1]:= (1/2)!
Out[1]= 1/2 Sqrt[Pi]
```

### Notes

`n!` and `Factorial[n]` compute exact integer factorials, promoting to GMP
bigints well before machine-word overflow, so `30!` is returned in full. The base
case `0! = 1` holds by convention. Half-integer arguments are evaluated through
the Gamma function, so `(1/2)! = Gamma[3/2] = Sqrt[Pi]/2`, printed as
`1/2 Sqrt[Pi]`. This connects the discrete factorial to its continuous Gamma
extension for non-integer inputs.
