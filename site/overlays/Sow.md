### Worked examples

```mathematica
In[1]:= Reap[Sow[1]; Sow[2]; Sow[3]]  (* {lastValue, {{everything sown}}} *)
```

```mathematica
In[1]:= Reap[Do[Sow[i^2], {i, 1, 4}]]  (* collect squares across a loop *)
```

```mathematica
In[1]:= Reap[Sow[1, "a"]; Sow[2, "b"]; Sow[3, "a"], _]  (* grouped by tag in first-encounter order *)
```

```mathematica
In[1]:= Sow[42]  (* returns its value; collection is a side effect *)
```

### Notes

`Sow[e]` records `e` into the nearest enclosing `Reap` and returns `e` unchanged — it is
transparent in an expression, so you can drop it into existing code to harvest
intermediate values. `Sow[e, tag]` routes `e` to the `Reap` whose pattern matches `tag`;
`Sow[e, {t1, t2, ...}]` records `e` once per tag.

`Reap[expr]` returns `{value, collected}` where `value` is the result of `expr` and
`collected` holds everything sown, grouped by tag in the order each tag was first seen.
Outside any `Reap`, `Sow` simply returns its argument (nothing is collected). Note that
`Sow` evaluates its argument normally — it is `Reap` that is `HoldFirst`. `Sow` is
`Protected`.
