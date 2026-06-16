### Worked examples

```mathematica
In[1]:= ReIm[3 + 4 I]
Out[1]= {3, 4}
```

It splits exact powers and quotients into their real/imaginary components —
`(2 + I)^3 = 2 + 11 I`, and a complex division reduces to integers:

```mathematica
In[1]:= ReIm[(2 + I)^3]
Out[1]= {2, 11}

In[2]:= ReIm[(3 + 4 I)/(1 - 2 I)]
Out[2]= {-1, 2}
```

On a transcendental argument it returns the numeric pair — here Euler's formula
`e^(iπ/4)` to 20 digits, the real and imaginary parts each `1/√2`:

```mathematica
In[1]:= ReIm[N[E^(I Pi/4), 20]]
Out[1]= {0.707106781186547524409, 0.707106781186547524395}
```

### Notes

`ReIm[z]` is shorthand for `{Re[z], Im[z]}`; a real-valued argument gives `{z, 0}`.
