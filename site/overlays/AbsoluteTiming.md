### Worked examples

```mathematica
In[1]:= AbsoluteTiming[Factorial[20]][[2]]  (* the result is the reproducible second element *)
```

```mathematica
In[1]:= Length[AbsoluteTiming[1 + 1]]  (* always a two-element {seconds, result} list *)
```

```mathematica
In[1]:= First[AbsoluteTiming[Pause[0]]] >= 0  (* the elapsed time is never negative *)
```

### Notes

`AbsoluteTiming[expr]` returns `{seconds, result}`, where `seconds` is the
elapsed **wall-clock** time measured from a monotonic clock. Only the second
element is reproducible — the timing varies between runs — so extract it with
`[[2]]`.

This, not `Timing`, is the right measurement for anything threaded. `Timing`
reports CPU time summed over threads, so the multithreaded reductions and
elementwise kernels, `Dot`, and the LAPACK-backed decompositions all read
roughly *core-count* times their true duration there; `AbsoluteTiming` reads the
time you actually waited. Because the clock is monotonic, a clock adjustment
mid-evaluation cannot produce a negative interval, and time spent sleeping in
`Pause` is included.
