### Worked examples

```mathematica
In[1]:= HoldComplete[1 + 1]  (* holds its argument, like Hold *)
```

```mathematica
In[1]:= HoldComplete[Sequence[1, 2], 3]  (* but unlike Hold it does not splice Sequence *)
```

```mathematica
In[1]:= Hold[Sequence[1, 2], 3]  (* for contrast, Hold does splice it *)
```

```mathematica
In[1]:= ReleaseHold[HoldComplete[1 + 1]]  (* ReleaseHold strips the wrapper and evaluates *)
```

### Notes

`HoldComplete[expr]` is the strongest hold wrapper: it carries the `HoldAllComplete`
attribute, so besides holding every argument unevaluated it also suppresses the upvalue
lookup and the `Sequence` / `Unevaluated`-stripping that an ordinary `Hold` still performs.
That is why `HoldComplete[Sequence[1, 2], 3]` keeps its `Sequence` intact where `Hold`
flattens it.

`ReleaseHold` removes the wrapper and lets the contents evaluate.
