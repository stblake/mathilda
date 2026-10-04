### Worked examples

```mathematica
In[1]:= Head[SessionTime[]]  (* always a Real *)
```

```mathematica
In[1]:= TrueQ[Positive[SessionTime[]]]  (* grows monotonically from kernel start-up *)
```

### Notes

`SessionTime[]` gives the wall-clock seconds elapsed since the current Mathilda
session began. It is measured from a monotonic clock captured at kernel
start-up — the same clock source as `AbsoluteTiming` — so it increases on every
call and includes time spent sleeping in `Pause`.

The exact value is non-deterministic (it depends on how long the session has
been running), so useful examples report its structure — its head, or its sign —
rather than a fixed number. Contrast `TimeUsed`, which counts only CPU time and
does not advance while the session is idle.
