### Worked examples

```mathematica
In[1]:= f[a, Sequence[b, c], d]  (* spliced into the enclosing argument list *)
```

```mathematica
In[1]:= {a, Sequence[b], c, Identity[d]}  (* a one-element Sequence is the identity *)
```

```mathematica
In[1]:= {a, Sequence[], c}  (* the empty Sequence evaporates *)
```

```mathematica
In[1]:= {a, b, g[x, y], h[w], g[z, y]} /. g -> Sequence
```

```mathematica
In[1]:= f[a, b, c] /. f[x__] -> x  (* BlankSequence binds to a Sequence object *)
```

### Notes

`Sequence[e1, e2, ...]` is a run of arguments that is automatically spliced into
the argument list of whatever function encloses it. Splicing happens structurally
during evaluation, **before** `Flat` / `Listable` / `Orderless`, so
`f[a, Sequence[b, c], d]` becomes `f[a, b, c, d]`. `Sequence[]` contributes
nothing and `Sequence[e]` contributes one argument.

A bare `Sequence[...]` with no enclosing function — including one stored in an
`OwnValue` — is left as a `Sequence` object and splices only when it reaches a call
site. It is the wrapper produced by `BlankSequence` / `BlankNullSequence` and by
`SlotSequence` (`##`). Splicing is suppressed under a head carrying `SequenceHold`
or `HoldAllComplete`.
