### Worked examples

```mathematica
In[1]:= MeanClusteringCoefficient[CompleteGraph[4]]  (* every neighbourhood is fully connected *)
```

```mathematica
In[1]:= MeanClusteringCoefficient[CycleGraph[5]]  (* a cycle of five has no triangles *)
```

```mathematica
In[1]:= MeanClusteringCoefficient[Graph[{1, 2, 3, 4}, {UndirectedEdge[1, 2], UndirectedEdge[2, 3], UndirectedEdge[1, 3], UndirectedEdge[3, 4]}]]  (* a triangle with a pendant vertex *)
```

```mathematica
In[1]:= MeanClusteringCoefficient[StarGraph[5]]  (* the leaves have degree 1, so their local coefficient is 0 *)
```

```mathematica
In[1]:= MeanClusteringCoefficient[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3], DirectedEdge[3, 1]}]]  (* a directed 3-cycle counts as a triangle *)
```

### Notes

This is the mean of `LocalClusteringCoefficient` over all vertices, where a vertex of degree `d` in `t` triangles scores `t / C(d, 2)`, or 0 when `d < 2`. The answer is exact, an integer or a rational number.

In the triangle-with-pendant example the local values are `1, 1, 1/3, 0`, whose mean is `7/12`. Mixed graphs, with both directed and undirected edges, are left unevaluated.
