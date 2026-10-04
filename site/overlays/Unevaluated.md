### Worked examples

```mathematica
In[1]:= Length[Unevaluated[Plus[5, 6, 7, 8]]]  (* Length sees the held Plus, not a sum *)
```

```mathematica
In[1]:= Length[Unevaluated[1 + 2 + 3]]
```

```mathematica
In[1]:= Length[Unevaluated[Sequence[a, b]]]  (* stripping happens after Sequence splicing *)
```

```mathematica
In[1]:= Hold[Unevaluated[1 + 2]]  (* a genuinely held slot keeps the wrapper *)
```

```mathematica
In[1]:= HoldComplete[Unevaluated[1 + 2]]
```

```mathematica
In[1]:= Attributes[Unevaluated]
```

### Notes

`f[Unevaluated[expr]]` makes an ordinary head `f` hold the single argument `expr`
that it would normally evaluate: the wrapper is stripped in the non-held slot, but
stripping it does **not** force evaluation of the exposed content for that step —
so `Length[Unevaluated[1 + 2 + 3]]` is `3`, because `Length` sees the held
`Plus[1, 2, 3]`.

The wrapper is **not** stripped in a slot that was already held — under `HoldAll`,
`HoldFirst`/`HoldRest`, or `HoldAllComplete` — because the argument was never going
to be evaluated: `Hold[Unevaluated[1 + 2]]` stays intact. Stripping happens after
`Sequence` flattening, so a `Sequence` directly inside `Unevaluated` survives into
the argument slot. `Unevaluated` itself carries `HoldAllComplete`, so as a
top-level expression it evaluates to itself.
