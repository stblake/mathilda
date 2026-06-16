### Worked examples

```mathematica
In[1]:= Chop[N[Sin[Pi]]]
Out[1]= 0
```

```mathematica
In[1]:= Chop[{1.5, 1.0*^-12, 3.2}]
Out[1]= {1.5, 0, 3.2}
```

```mathematica
In[1]:= Chop[N[Exp[I Pi] + 1]]
Out[1]= 0
```

```mathematica
In[1]:= Chop[3.0 + 1.0*^-15 I]
Out[1]= 3.0
```

```mathematica
In[1]:= Chop[1.0*^-15 + 2.5 I]
Out[1]= 0.0 + 2.5*I
```

### Notes

`Chop[expr]` replaces approximate reals within `delta` of zero (default
`10^-10`) by the exact integer `0`, walking the whole expression tree. It is the
standard tool for cleaning numerical noise: `Chop[N[Sin[Pi]]]` and Euler's
identity `Chop[N[Exp[I Pi] + 1]]` both collapse to a clean `0`. For machine
complex numbers it chops parts independently — dropping a negligible imaginary
part recovers a pure real (`3.0`), while a negligible real part leaves a
`Complex[0., im]` that preserves the machine-complex shape. Exact numbers and
symbolic constants pass through untouched.
