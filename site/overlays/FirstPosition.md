### Worked examples

```mathematica
In[1]:= FirstPosition[{a, b, c, b}, b]
```

```mathematica
In[1]:= FirstPosition[{{1, 2}, {3, 4}}, 4]  (* a nested position, depth first *)
```

```mathematica
In[1]:= FirstPosition[{1, 2, 3}, 5]
```

```mathematica
In[1]:= FirstPosition[{1, 2, 3}, 5, None]  (* the held default is returned when nothing matches *)
```

```mathematica
In[1]:= FirstPosition[<|"a" -> 1, "b" -> 2|>, 2]  (* over an association the position is a key *)
```

```mathematica
In[1]:= FirstPosition[{1, 2, 3}, _Integer, None, {1}]  (* restricted to level 1 *)
```

### Notes

`FirstPosition[expr, pattern]` gives the position — a list of indices — of the
first subexpression matching `pattern` in depth-first order, or
`Missing["NotFound"]` if there is none. `FirstPosition[expr, pattern, default]`
returns `default` instead; `default` is held and evaluated only when it is
actually returned. A fourth argument is a level specification.

It delegates to `Position` (with the match count fixed at 1), so it inherits
`Position`'s traversal, default level spec `{0, Infinity}`, `Heads -> True`
handling, and association value → `Key[...]` remapping exactly. Because the
default spec includes level 0, a pattern that matches the whole of `expr` gives
the position `{}`.
