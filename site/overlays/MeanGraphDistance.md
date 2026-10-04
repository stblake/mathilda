### Worked examples

```mathematica
In[1]:= MeanGraphDistance[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5}]]  (* exact for an unweighted graph *)
```

```mathematica
In[1]:= MeanGraphDistance[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight -> {1, 4}]]  (* weighted: a machine real *)
```

### Notes

The value is the mean distance over all ordered pairs of distinct vertices. It is exact (an
`Integer` or `Rational`) for an unweighted graph and a machine `Real` for a weighted one.

It is `Infinity` unless the graph is connected (strongly connected, for a directed graph), and
is undefined — left unevaluated — for a graph with a single vertex.
