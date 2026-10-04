### Worked examples

```mathematica
In[1]:= HoldForm[1 + 1]  (* held, but the wrapper prints invisibly *)
```

```mathematica
In[1]:= FullForm[HoldForm[1 + 1]]  (* the wrapper is really there *)
```

```mathematica
In[1]:= ReleaseHold[HoldForm[1 + 1]]  (* ReleaseHold strips it and evaluates *)
```

### Notes

`HoldForm[expr]` holds `expr` unevaluated exactly like `Hold`, but the printer renders it
as just `expr` — the wrapper is invisible in output. It is the tool for displaying an
expression in unevaluated form while keeping it a genuine held expression, as `FullForm`
reveals.

`HoldForm` carries `HoldAll` and is Protected; `ReleaseHold` removes it.
