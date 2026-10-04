### Worked examples

```mathematica
In[1]:= I^2  (* the defining relation: the square of the imaginary unit is -1 *)
```

```mathematica
In[1]:= Sqrt[-1]  (* the principal square root of -1 folds back to I *)
```

```mathematica
In[1]:= Re[I^2]  (* the real part of the square *)
```

```mathematica
In[1]:= Exp[I Pi]  (* Euler's identity *)
```

```mathematica
In[1]:= Exp[I Pi/2]  (* a quarter turn around the unit circle *)
```

```mathematica
In[1]:= (1 + I) (1 - I)  (* z times its conjugate is the squared modulus *)
```

```mathematica
In[1]:= Conjugate[I]  (* conjugation flips the sign of the imaginary part *)
```

```mathematica
In[1]:= {Abs[I], Arg[I]}  (* modulus one, argument a right angle *)
```

```mathematica
In[1]:= ComplexExpand[(1 + I)^2]  (* expanding a complex power *)
```

```mathematica
In[1]:= FullForm[I]  (* internally the imaginary unit is a Complex atom *)
```

### Notes

`I` is the imaginary unit √(−1). It is not stored as a mathematical constant in its
own right: the symbol `I` carries an OwnValue rewriting it to `Complex[0, 1]`, so
`FullForm[I]` is `Complex[0, 1]` and `Head[I]` is `Complex`. Its only attribute is
`Protected` (it is *not* `Constant`), yet `NumericQ[I]` is `True` and `D[I, x]` is `0`
because it resolves to a number. Every complex identity above — `I^2 == -1`,
`Sqrt[-1] == I`, `Conjugate[I] == -I`, `Abs[I] == 1`, `Arg[I] == Pi/2`, and the Euler
relations `Exp[I Pi] == -1`, `Exp[I Pi/2] == I` — falls out of the generic complex
arithmetic. Since the real and imaginary parts are the exact integers `0` and `1`,
`N[I]` is simply `0. + 1. I`: there is no arbitrary-precision form of `I` to request.
