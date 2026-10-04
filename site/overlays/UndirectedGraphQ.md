### Worked examples

```mathematica
In[1]:= UndirectedGraphQ[CycleGraph[5]]  (* every edge of a cycle graph is undirected *)
```

```mathematica
In[1]:= UndirectedGraphQ[Graph[{1 -> 2}]]  (* a single directed edge *)
```

```mathematica
In[1]:= UndirectedGraphQ[Graph[{1 <-> 2, 2 <-> 3}]]  (* explicitly undirected edges *)
```

```mathematica
In[1]:= UndirectedGraphQ[CompleteGraph[4]]  (* the generators build undirected graphs *)
```

```mathematica
In[1]:= UndirectedGraphQ[7]  (* a non-graph argument is False *)
```

### Notes

`UndirectedGraphQ[g]` is `True` precisely when `g` has no directed edges, which includes
the edgeless case: a graph with no edges is undirected because nothing in it is oriented.
A purely directed graph, or one that mixes the two edge kinds, gives `False`.

The test is a single memoized directed-edge count, so it is `O(1)` after the first query
on a given graph. A non-graph argument yields `False` rather than remaining unevaluated.
