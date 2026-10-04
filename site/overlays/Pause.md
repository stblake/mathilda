### Worked examples

```mathematica
In[1]:= Pause[0] === Null  (* Pause returns Null; zero waits not at all *)
```

```mathematica
In[1]:= Pause[1/100] === Null  (* a rational (or any NumericQ) duration is accepted *)
```

### Notes

`Pause[n]` blocks for at least `n` seconds of wall-clock time, then returns
`Null`. The argument may be any non-negative number — integer, real, rational
(`Pause[1/4]`), or a `NumericQ` symbolic form (`Pause[Pi]`); zero or negative
returns immediately, and a non-numeric argument leaves `Pause[x]` unevaluated.

The wait is a `nanosleep`, which consumes no CPU. So the elapsed time is counted
by `AbsoluteTiming` and `SessionTime` but **not** by `Timing` or `TimeUsed`:
`Timing[Pause[1]]` reports about `0`, while `AbsoluteTiming[Pause[1]]` reports
about `1`. The sleep resumes across signal interruptions, so the full duration is
always observed. Accuracy is bounded by the timer granularity `$TimeUnit`.
