### Worked examples

```mathematica
In[1]:= BesselJ[0, 0]
Out[1]= 1

In[2]:= BesselJ[1, 0]
Out[2]= 0
```

Half-integer orders close in elementary functions:

```mathematica
In[1]:= BesselJ[1/2, z]
Out[1]= Sin[z] Sqrt[2/(Pi z)]
```

The Frobenius series at the origin:

```mathematica
In[1]:= Series[BesselJ[0, x], {x, 0, 6}]
Out[1]= 1 - 1/4 x^2 + 1/64 x^4 - 1/2304 x^6 + O[x]^7
```

High-precision and complex evaluation; the first input is `J_0` at its first zero, which numerically returns essentially zero:

```mathematica
In[1]:= N[BesselJ[0, 1], 40]
Out[1]= 0.7651976865579665514497175261026632209093

In[2]:= N[BesselJ[0, 10 + 5 I], 30]
Out[2]= -17.78959112945037151834426180967 + 0.2007116167212048509818027697064*I
```

### Notes

`BesselJ[n, z]` is the Bessel function of the first kind, regular at the origin, with `J_0(0) = 1` and `J_n(0) = 0` for integer `n != 0`. Real and complex order and argument evaluate at machine or MPFR precision; `D[BesselJ[n, z], z] = (BesselJ[n-1, z] - BesselJ[n+1, z])/2`. There is a branch cut along the negative real axis for non-integer order. Listable.
