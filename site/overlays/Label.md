### Worked examples

```mathematica
In[1]:= Module[{k = 1}, Label[a]; k = 2 k; If[k < 16, Goto[a]]; k]  (* Label marks the target that Goto returns to *)
```

```mathematica
In[1]:= Label[done]  (* as a bare statement a Label evaluates to Null *)
```

### Notes

`Label[tag]` marks a point that `Goto[tag]` can jump to. It must appear as an
explicit element of a `CompoundExpression` — it is the literal `Label[tag]` node
that the compound expression scans for when it consumes a `Goto` sentinel.

Evaluated in ordinary top-to-bottom flow a `Label` is a no-op worth `Null`, so
reaching one by falling through simply continues to the next statement; only a
`Goto` gives it an effect.
