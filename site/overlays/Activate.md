### Worked examples

```mathematica
In[1]:= Activate[Inactive[Plus][2, 3]]  (* reactivates the inert head, so the sum is taken *)
```

```mathematica
In[1]:= Activate[Inactive[Integrate][x, x]]  (* an inert integral becomes a real one and evaluates *)
```

### Notes

`Activate[expr]` is the inverse of `Inactive`: it walks the expression replacing every
inert head `Inactive[h]` by `h` and then lets the evaluator re-run, so the reactivated
heads finally fire. It is the counterpart to holding a computation inert with `Inactive`.

The rewrite is recursive and reaches every nesting level, so a deeply inert tree is
reactivated in a single pass.
