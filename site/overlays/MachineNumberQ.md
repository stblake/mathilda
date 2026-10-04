### Worked examples

```mathematica
In[1]:= MachineNumberQ[Sin[1000.]]  (* a finite machine double *)
```

```mathematica
In[1]:= MachineNumberQ[Exp[1000.]]  (* overflows to +inf, so not a machine number *)
```

```mathematica
In[1]:= MachineNumberQ[N[Pi, 30]]  (* arbitrary-precision MPFR, not machine *)
```

```mathematica
In[1]:= MachineNumberQ[1.0 + 2.0 I]  (* two finite machine reals *)
```

```mathematica
In[1]:= MachineNumberQ[1 + 2 I]  (* exact Gaussian integer, not machine *)
```

### Notes

`MachineNumberQ[expr]` is `True` for a machine-precision (IEEE double) real, or a
`Complex` whose real and imaginary parts are both finite machine reals. It draws
three distinctions that `NumberQ` does not: an exact number (integer, bigint,
rational) is **not** a machine number; an arbitrary-precision MPFR real such as
`N[Pi, 30]` is **not**; and a non-finite double (an overflow to `inf`, or a `NaN`)
is **not**, so a computation that silently overflowed can be caught with this
test.
