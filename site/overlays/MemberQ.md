### Worked examples

```mathematica
In[1]:= MemberQ[{1, 2, 3}, 2]
Out[1]= True
```

```mathematica
In[1]:= MemberQ[{1, 2, 3}, _Integer]
Out[1]= True
```

```mathematica
In[1]:= MemberQ[{x^2, y^3, z}, _^_]
Out[1]= True
```

```mathematica
In[1]:= MemberQ[{{1, 2}, {3, 4}}, 3, {2}]
Out[1]= True
```

```mathematica
In[1]:= MemberQ[#, 0] & /@ {{1, 2}, {0, 3}}
Out[1]= {False, True}
```

### Notes

`MemberQ[list, form]` tests whether any element of `list` matches `form`. The
second argument is a full pattern, not just a literal value — `MemberQ[list,
_Integer]` checks for any integer and `MemberQ[{x^2, y^3, z}, _^_]` finds the
first element stored as a `Power`. The optional level specification controls how
deep to look: `MemberQ[{{1, 2}, {3, 4}}, 3, {2}]` searches at level 2 (inside the
sub-lists) and so detects the `3`, whereas the default level 1 would not. Because
`MemberQ` is an ordinary predicate, it maps cleanly over collections, as in the
last example which flags which rows contain a zero.
