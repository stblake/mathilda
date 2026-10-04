### Worked examples

```mathematica
In[1]:= RationalQ[3/4]  (* an exact Rational[3, 4] *)
```

```mathematica
In[1]:= RationalQ[7]  (* an integer is rational *)
```

```mathematica
In[1]:= RationalQ[2.5]  (* a machine real is not, even with a rational value *)
```

```mathematica
In[1]:= RationalQ[Pi]  (* an irrational constant is not rational *)
```

### Notes

`RationalQ[expr]` is `True` for exactly the exact rationals: any integer (`Integer`
or big `BigInt`) and any `Rational[p, q]`. It is a representation test, so reals and
MPFR numbers are `False` regardless of their value, and irrational constants and
symbolic expressions are `False` as well. Non-single-argument calls are left
unevaluated.
