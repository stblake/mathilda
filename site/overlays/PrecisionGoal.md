### Worked examples

```mathematica
In[1]:= Head[PrecisionGoal]  (* a bare option symbol, not a function *)
```

```mathematica
In[1]:= FullForm[PrecisionGoal -> 10]  (* it only sits on the left of a rule *)
```

```mathematica
In[1]:= NIntegrate[1/x, {x, 1, 2}, PrecisionGoal -> 10]  (* consumed by the numerical op *)
```

### Notes

`PrecisionGoal` is the option shared by the numerical operations (`NIntegrate`, `NDSolve`,
`NSum`, …) that says how many digits of **relative** precision to seek. It is an inert
option-name symbol — no builtin, no DownValues — so it does nothing on its own; each
numerical subsystem reads it and turns it into a tolerance.

With `PrecisionGoal -> p` and `AccuracyGoal -> a` the operation seeks an error below
`10^-a + |x| 10^-p` in a value of size `x`, so `PrecisionGoal` bounds the relative term.
The value may be a digit count, `Automatic` (the default, two digits below
`WorkingPrecision`), or `Infinity` (disable the relative criterion). If the goal cannot be
met at the resource cap, a `Head::accgl` message is emitted and the best approximation
returned; keep `WorkingPrecision` at least as large as the goal.
