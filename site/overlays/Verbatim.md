### Worked examples

```mathematica
In[1]:= {x_, x} /. Verbatim[x_] -> matched  (* only the literal x_ is rewritten, not the symbol x *)
```

```mathematica
In[1]:= MatchQ[x_, Verbatim[x_]]  (* the literal pattern expression matches itself *)
```

```mathematica
In[1]:= MatchQ[5, Verbatim[x_]]  (* but an ordinary value does not *)
```

```mathematica
In[1]:= Count[{a + b, x_ + y_, 1 + 2}, Verbatim[x_ + y_]]  (* find a literal pattern sitting inside data *)
```

### Notes

`Verbatim[expr]` matches `expr` taken literally — the pattern constructs inside it
(`Blank`, `Pattern`, ...) are not interpreted. So `Verbatim[x_]` matches only the
expression `x_` itself, which is what lets you search for, or rewrite, a pattern
*as data*. The matcher implements it as a structural-equality test (`expr_eq`), so
it binds no variables; it is transparent in the same way as `HoldPattern`, but
where `HoldPattern` keeps its argument interpretable as a pattern, `Verbatim`
freezes it to a literal.
