### Worked examples

```mathematica
In[1]:= Inactive[Plus][2, 3]  (* the head is held inert, so the sum is not taken *)
```

```mathematica
In[1]:= Inactive[Plus][1 + 1, 3]  (* but the arguments still evaluate *)
```

```mathematica
In[1]:= Inactive[Integrate][x, x]  (* a symbolic, inert integral *)
```

```mathematica
In[1]:= Activate[Inactive[Integrate][x, x]]  (* Activate turns it back on *)
```

### Notes

`Inactive[h]` wraps a head so that its own rules are suppressed: `Inactive[h][args]` has
the compound head `Inactive[h]`, which has no rules of its own, so the call stays a
symbolic fixed point. It is how a computation like an integral is kept inert for display
or later manipulation.

Only the head is held — the arguments still evaluate, so `Inactive[Plus][1 + 1, 3]` first
reduces to `Inactive[Plus][2, 3]`. `Activate` reverses the wrapper, and `D` already knows
how to differentiate an inert `Inactive[Integrate]`.
