### Worked examples

```mathematica
In[1]:= Module[{i = 0, s = 0}, Label[top]; i++; s += i; If[i < 5, Goto[top]]; s]  (* a backward jump forms a loop *)
```

```mathematica
In[1]:= Module[{k = 1}, Label[a]; k = 2 k; If[k < 16, Goto[a]]; k]  (* doubling until the guard fails *)
```

```mathematica
In[1]:= Goto[nowhere]  (* with no matching Label it is left inert *)
```

### Notes

`Goto[tag]` transfers control to the `Label[tag]` in the `CompoundExpression` the
`Goto` appears in directly, then in enclosing ones — a forward jump skips the
statements in between, a backward jump forms a loop. Like `Throw`, it propagates
by a sentinel through the evaluator's normal return path, so a `Goto` fired inside
a nested call (an `If` branch, as above) still reaches the enclosing
`CompoundExpression`.

A `Goto` loop has no artificial iteration cap; termination is the program's
responsibility, exactly as for `While`. A `Goto[tag]` that reaches top level with
no matching `Label` anywhere prints `Goto::nolabel` and returns the inert
`Goto[tag]` node.
