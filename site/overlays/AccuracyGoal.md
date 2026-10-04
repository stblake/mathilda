### Worked examples

```mathematica
In[1]:= Head[AccuracyGoal]  (* a bare option symbol, not a function *)
```

```mathematica
In[1]:= FullForm[AccuracyGoal -> 8]  (* it only sits on the left of a rule *)
```

```mathematica
In[1]:= NIntegrate[Exp[-x^2], {x, 0, 1}, AccuracyGoal -> 8]  (* consumed by the numerical op *)
```

### Notes

`AccuracyGoal` is the option shared by the numerical operations (`NIntegrate`, `NDSolve`,
`NSum`, `FindRoot`, `NSolve`, …) that says how many digits of **absolute** accuracy to
seek. It is an inert option-name symbol — no builtin, no DownValues — so it does nothing on
its own; each numerical subsystem reads it and turns it into a tolerance.

With `AccuracyGoal -> a` and `PrecisionGoal -> p` the operation seeks an error below
`10^-a + |x| 10^-p` in a value of size `x`, so `AccuracyGoal` bounds the absolute term.
The value may be a digit count, `MachinePrecision` (the default), `Automatic`, or
`Infinity` (disable the absolute criterion). If the goal cannot be met at the resource cap
the operation emits a `Head::accgl` message and returns its best approximation; keep
`WorkingPrecision` at least as large as the goal.
