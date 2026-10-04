### Worked examples

```mathematica
In[1]:= Trace[1 + 1]
```

```mathematica
In[1]:= Trace[5]  (* an inert atom needs no rewriting *)
```

```mathematica
In[1]:= Trace[2^3 + 4^2 + 1]  (* nested, mirroring the evaluation structure *)
```

```mathematica
In[1]:= Trace[1 + 2 + 3, _Integer]  (* filtered to steps matching a pattern *)
```

```mathematica
In[1]:= Trace[Nest[f, x, 3], _f]
```

### Notes

`Trace[expr]` returns a **nested** list of the forms `expr` passes through while
evaluating, mirroring the evaluator's own recursion: each argument sub-evaluation
that takes a step becomes a sublist, and the reassembled intermediate form appears
as a step. An expression that needs no rewriting (an inert atom, a normal form)
traces to `{}`. A builtin's internal work and `Listable` threading are shown as a
single atomic step, matching Mathematica — `Range[10]` is one step, not ten.

`Trace[expr, form]` keeps only the step leaves whose expression matches the
pattern `form`, flattening the nesting into a plain list. `form` is held, so
pattern literals such as `_Integer` or `f[_]` may be written directly. Each step
is returned wrapped in `HoldForm`, so the result prints transparently yet stays
inert and does not re-evaluate.
