### Worked examples

```mathematica
In[1]:= GraphRadius[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5}]]  (* the central vertex has eccentricity 2 *)
```

```mathematica
In[1]:= gw = Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight -> {1, 4}];
```

```mathematica
In[1]:= GraphRadius[gw]  (* weighted: a machine real *)
```

### Notes

The radius is the minimum eccentricity over all vertices (the eccentricity of a vertex being its
greatest distance to any other). It is an exact `Integer` for an unweighted graph and a machine
`Real` for a weighted one.

It is `Infinity` unless the graph is connected (strongly connected, for a directed graph).
`GraphCenter` returns the vertices that attain this minimum eccentricity.
