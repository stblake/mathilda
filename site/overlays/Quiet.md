### Worked examples

```mathematica
In[1]:= Quiet[1/0]  (* the Power::infy message is silenced; the value is unchanged *)
```

```mathematica
In[1]:= Quiet[Log[0]]  (* the value comes back with no printed warning *)
```

### Notes

`Quiet[expr]` evaluates `expr` and returns its value with any messages suppressed.
It is `HoldAll`, so the suppression is in force *during* the evaluation. The
optional second argument (a message name or list) is accepted and ignored — all
messages are suppressed, which is a harmless superset of any requested set.

Suppression silences only the *printing*. A message still fires while quieted, so
an enclosing `Check` still sees it — which is exactly what the common
`Quiet[Check[expr, failexpr]]` pattern relies on.
