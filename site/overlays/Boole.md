### Worked examples

```mathematica
In[1]:= Boole[3 > 2]
Out[1]= 1
```

```mathematica
In[1]:= Boole[{True, False, True}]
Out[1]= {1, 0, 1}
```

```mathematica
In[1]:= Table[Boole[PrimeQ[n]], {n, 1, 12}]
Out[1]= {0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0}
```

```mathematica
In[1]:= Sum[Boole[GCD[k, 10] == 1], {k, 1, 10}]
Out[1]= 4
```

```mathematica
In[1]:= Sum[Boole[Mod[k, 3] == 0] k^2, {k, 1, 10}]
Out[1]= 126
```

### Notes

`Boole[expr]` is the Iverson bracket: it yields `1` when `expr` is `True` and
`0` when it is `False`, and stays unevaluated otherwise. Being `Listable`, it
threads element-wise over lists. Combined with `Sum` it turns logical predicates
into counting and conditional-summation devices: `Sum[Boole[GCD[k, 10] == 1], {k, 1, 10}]`
counts the integers up to 10 coprime to 10 (Euler's totient phi(10) = 4), and
weighting the bracket by `k^2` sums squares restricted to a residue class.
