### Worked examples

```mathematica
In[1]:= MaximalBy[{1, -5, 3, -2}, Abs]  (* the element with the largest absolute value *)
```

```mathematica
In[1]:= MaximalBy[{"a", "bbb", "cc", "ddd"}, StringLength]  (* all ties, in input order *)
```

```mathematica
In[1]:= MaximalBy[Abs][{1, -5, 3}]  (* operator form MaximalBy[f][list] *)
```

### Notes

`MaximalBy[list, f]` returns the element (or elements) of `list` for which `f` is
largest, measured by Mathilda's canonical order. All elements tying for the
maximum are returned together, in their original order, so the result is always a
list — here `{"bbb", "ddd"}` both have length 3. It selects an extreme rather than
sorting. The operator form `MaximalBy[f]` applies to a collection supplied later;
over an association the entries whose *value* maximises `f` are returned.
