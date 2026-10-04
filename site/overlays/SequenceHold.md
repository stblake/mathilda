### Worked examples

```mathematica
In[1]:= Hold[Sequence[a, b]]  (* HoldAll alone does not stop the splice *)
```

```mathematica
In[1]:= SetAttributes[h, {HoldAll, SequenceHold}]; h[a, Sequence[b, c]]  (* now it is held *)
```

```mathematica
In[1]:= MemberQ[Attributes[Set], SequenceHold]  (* assignment heads carry it *)
```

```mathematica
In[1]:= splice[x_] := Sequence[x, x, x]; {a, splice[b], c}  (* splices at the call site *)
```

### Notes

`SequenceHold` is an **attribute**, not a function: it tells the evaluator not to
flatten `Sequence[...]` objects that appear among a function's arguments.
`HoldAll` on its own does *not* prevent the splice — `Hold[Sequence[a, b]]` gives
`Hold[a, b]` — so `SequenceHold` is the specific lever that stops it, and
`HoldAllComplete` implies it.

The assignment and rule heads `Set`, `SetDelayed`, `Rule`, and `RuleDelayed` all
carry `SequenceHold`. That is what lets a definition store a `Sequence` on its
right-hand side and have it splice only when the defined symbol is later used — the
`splice[x_] := Sequence[x, x, x]` idiom. Any user head that sets the attribute is
honoured the same way.
