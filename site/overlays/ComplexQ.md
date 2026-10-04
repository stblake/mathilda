### Worked examples

```mathematica
In[1]:= ComplexQ[2 + 3 I]  (* stored as Complex[2, 3], so the head test passes *)
```

```mathematica
In[1]:= ComplexQ[2.0 + 3.0 I]  (* inexact parts are still a Complex *)
```

```mathematica
In[1]:= ComplexQ[2.0]  (* a bare real is not Complex *)
```

```mathematica
In[1]:= ComplexQ[5]  (* nor is an integer *)
```

### Notes

`ComplexQ[expr]` tests one thing: whether `expr` has head `Complex`. It is a
*syntactic* test on the representation, not a check of mathematical type — a real
number is stored with a real head, so `ComplexQ` is `False` for it regardless of
value. The parts may be exact or inexact; only the head matters. Anything that is
not a single-argument call is left unevaluated.
