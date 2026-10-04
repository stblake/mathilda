### Worked examples

```mathematica
In[1]:= FirstCase[{1, "a", 2, "b"}, _String]  (* the first element matching the pattern *)
```

```mathematica
In[1]:= FirstCase[{1, 2, 3, 4}, x_ /; x > 2]  (* a conditional pattern *)
```

```mathematica
In[1]:= FirstCase[{1, 2, 3, 4}, x_?EvenQ -> x^2]  (* a transformation rule returns the rewritten match *)
```

```mathematica
In[1]:= FirstCase[{1, 2, 3}, _String, None]  (* no match, so the supplied default *)
```

```mathematica
In[1]:= FirstCase[{1, 2, 3}, _String]  (* no match and no default: Missing["NotFound"] *)
```

### Notes

`FirstCase[expr, pattern]` gives the first element of `expr` matching `pattern`,
the single-element companion to `Cases`. When `pattern` is a transformation rule
`patt -> rhs`, the result is the rewritten first match rather than the raw element
(here `4` is `2^2`). With no match it returns `Missing["NotFound"]`, or the
`default` passed as a third argument.
