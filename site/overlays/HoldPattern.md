### Worked examples

```mathematica
In[1]:= HoldPattern[1 + 1]
Out[1]= HoldPattern[1 + 1]
```

The usual reason to reach for `HoldPattern` is on the left-hand side of a
rule: it keeps a structural pattern from being evaluated away before matching.
Here `HoldPattern[p_ + q_]` matches each unevaluated symbolic sum and rewrites
it as a product:

```mathematica
In[1]:= {a + b, c + d} /. HoldPattern[p_ + q_] -> p*q
Out[1]= {a b, c d}
```

It lets a rule target a still-unevaluated head, even one the evaluator would
normally leave inert, such as a symbolic `Integrate`:

```mathematica
In[1]:= Integrate[f[x], x] /. HoldPattern[Integrate[a_, b_]] -> done
Out[1]= done
```

`ReleaseHold` strips one wrapping layer, so the held expression finally
evaluates:

```mathematica
In[1]:= ReleaseHold[HoldPattern[2 + 3]]
Out[1]= 5
```

### Notes

`HoldPattern[expr]` is equivalent to `expr` for matching purposes but keeps `expr` unevaluated, so the LHS of a rule or assignment is not simplified before it is used as a pattern. It has attributes `{HoldAll, Protected}` and is removed by one layer of `ReleaseHold`.
