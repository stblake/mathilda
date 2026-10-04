### Worked examples

```mathematica
In[1]:= Head[TimeUsed[]]  (* always a Real *)
```

```mathematica
In[1]:= TimeUsed[] >= 0  (* CPU seconds consumed so far this session *)
```

### Notes

`TimeUsed[]` gives the total CPU seconds consumed so far in the current session,
read from the same processor clock `Timing` brackets. Like `Timing`, it does
**not** advance during `Pause` or other idle waits — only work that actually
burns CPU moves it — and it counts CPU time summed over threads.

The value is non-deterministic (it grows with the work the session has done), so
examples show its structure rather than a fixed number. For elapsed real time
since the session began, including idle and `Pause` time, use `SessionTime`
instead; its resolution is `$TimeUnit`.
