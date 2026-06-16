### Worked examples

```mathematica
In[1]:= Sort[{3, 1, 2}]
Out[1]= {1, 2, 3}

In[2]:= Sort[{5, 3, 8, 1}, Greater]
Out[2]= {8, 5, 3, 1}
```

```mathematica
In[1]:= Sort[{x^2, x, 1, x^3}]
Out[1]= {1, x, x^2, x^3}

In[2]:= Sort[{"banana", "apple", "cherry"}]
Out[2]= {"apple", "banana", "cherry"}
```

```mathematica
In[1]:= Sort[{{2, 1}, {1, 3}, {1, 2}}]
Out[1]= {{1, 2}, {1, 3}, {2, 1}}

In[2]:= Sort[Range[10], (Mod[#1, 3] < Mod[#2, 3]) &]
Out[2]= {3, 9, 6, 10, 1, 7, 4, 2, 8, 5}
```

### Notes

`Sort[list]` orders elements by Mathilda's canonical ordering, which compares numbers numerically, strings lexicographically, and structured expressions component-by-component (so the nested lists sort by first element, then second). `Sort[list, p]` uses an ordering predicate `p[a, b]` instead: `Greater` reverses to descending order, and a pure function such as `Mod[#1, 3] < Mod[#2, 3]` groups by residue class. The sort is stable, so equal-ranked elements keep their original relative order.
