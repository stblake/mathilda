### Worked examples

```mathematica
In[1]:= Factorial2[7]
Out[1]= 105

In[2]:= Factorial2[8]
Out[2]= 384

In[3]:= 0!!
Out[3]= 1
```

```mathematica
In[1]:= Factorial2[19] Factorial2[20]
Out[1]= 2432902008176640000

In[2]:= Factorial2[19] Factorial2[20] == 20!
Out[2]= True
```

### Notes

`Factorial2[n]` (also typeset `n!!`) is the double factorial: it multiplies down in steps of 2, so `7!! = 7*5*3*1 = 105` and `8!! = 8*6*4*2 = 384`. By convention `0!! = (-1)!! = 1`. The interleaving identity `(n-1)!! · n!! = n!` is exact here: the product of the odd double factorial `19!!` and the even double factorial `20!!` reproduces `20!` term for term.
