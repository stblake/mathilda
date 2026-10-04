### Worked examples

```mathematica
In[1]:= Implies[True, q]  (* a true premise reduces to the conclusion *)
```

```mathematica
In[1]:= Implies[False, q]  (* a false premise implies anything *)
```

```mathematica
In[1]:= Implies[p, False]  (* p implies False is the negation of p *)
```

```mathematica
In[1]:= Implies[p, p]  (* a statement implies itself *)
```

```mathematica
In[1]:= LogicalExpand[Implies[p, q]]  (* the definition, !p || q *)
```

### Notes

`Implies[p, q]` is material implication, `p ⟹ q`, logically `!p || q`. It takes
exactly two arguments and is neither associative nor commutative, so (unlike `And`
/`Or`/`Xor`/`Equivalent`) it is not `Flat` or `Orderless`.

It simplifies the literal and structural cases — `Implies[False, q]` and
`Implies[p, True]` are `True`, `Implies[True, q]` is `q`, `Implies[p, False]` is
`Not[p]`, and `Implies[p, p]` is `True` — and otherwise stays symbolic.
`LogicalExpand` and `Reduce` rewrite it to `!p || q` when a boolean normal form is
wanted.
