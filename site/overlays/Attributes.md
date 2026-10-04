### Worked examples

```mathematica
In[1]:= Attributes[Plus]  (* the full bundle an arithmetic head carries *)
```

```mathematica
In[1]:= Attributes[Sin]  (* a listable numeric function *)
```

```mathematica
In[1]:= Attributes[Hold]  (* HoldAll collapses the HoldFirst and HoldRest bits into one token *)
```

### Notes

`Attributes[sym]` returns the sorted list of attribute symbols set on `sym`, read from its
attribute bitflags. The `HoldFirst | HoldRest` pair is reported as the single token
`HoldAll`; `Protected` marks every builtin.

`Attributes` holds its argument (`HoldAll`), so the symbol is not evaluated first — you
get the attributes of `sym` itself, not of its value. Set and clear them with
`SetAttributes` and `ClearAttributes`.
