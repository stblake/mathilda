### Worked examples

```mathematica
In[1]:= ClearAll[g]; SetAttributes[g, HoldRest]; Attributes[g]  (* set the attribute on a head *)
```

```mathematica
In[1]:= ClearAll[g]; SetAttributes[g, HoldRest]; g[1 + 1, 2 + 2]  (* first arg evaluates, the rest are held *)
```

```mathematica
In[1]:= MemberQ[Attributes[If], HoldRest]  (* If holds its branches *)
```

### Notes

`HoldRest` is an attribute: when a head has it, every argument *except the first* is kept
unevaluated. In `g[1 + 1, 2 + 2]` above, `1 + 1` evaluates to `2` while `2 + 2` is left
as written, giving `g[2, 2 + 2]`.

It is the complement of `HoldFirst` and the partial form of `HoldAll`. Control
constructs such as `If` carry it so a branch does not evaluate until the condition picks
it; `RuleDelayed` carries it too. Use `Evaluate` to force a held argument to evaluate in
a controlled way. Like all attributes it is global to the head and is cleared by
`ClearAll` or `ClearAttributes`.
