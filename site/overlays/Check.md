### Worked examples

```mathematica
In[1]:= Check[2 + 2, failed]  (* no message fires, so the value passes through *)
```

```mathematica
In[1]:= Check[1/0, failed]  (* the division emits Power::infy, so failexpr is returned *)
```

```mathematica
In[1]:= Quiet[Check[1/0, caughtit]]  (* detect the failure without printing its message *)
```

### Notes

`Check[expr, failexpr]` returns `failexpr` if *any* message is generated while
`expr` is evaluated, otherwise the value of `expr`. It is `HoldAll`, so `expr` is
evaluated under `Check`'s watch rather than before it.

`Check` keys off whether a message *fired*, which is independent of whether it was
*printed* — so `Quiet[Check[expr, failexpr]]` is the idiom for detecting a failure
silently. A `Throw` inside `expr` is not a message: it propagates straight out of
`Check` rather than tripping the failure branch.
