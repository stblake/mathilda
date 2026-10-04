### Worked examples

```mathematica
In[1]:= Log10[1000]  (* powers of ten come back exact, via Log[10, z]'s power detection *)
```

```mathematica
In[1]:= Log10[{1, 10, 100, 1000}]  (* Listable, so a whole list maps at once *)
```

```mathematica
In[1]:= N[Log10[2], 30]  (* an arbitrary-precision value routes through the MPFR logarithm *)
```

```mathematica
In[1]:= Log10[-1.]  (* a negative real lands off the real axis, exactly as Log does *)
```

```mathematica
In[1]:= D[Log10[x], x]  (* differentiated through the Log[10, z] = Log[z]/Log[10] definition *)
```

### Notes

`Log10[z]` is `Log[10, z] = Log[z]/Log[10]`, so every rule `Log` knows — exact
powers, the branch cut on the negative real axis, arbitrary precision — is
inherited. The only thing `Log10` adds is the dedicated libm `log10` on a
positive machine real, which keeps `Log10[1000.]` exactly `3.` where the quotient
of two rounded logarithms would land an ulp away.

The real-valued `NDArray` kernel maps a packed buffer element-wise and compiles
at scalar and rank-1 shapes, so `Log10` in a `Compile[]`d body or over a packed
column runs on the buffer rather than one boxed number at a time.
