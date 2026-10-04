### Worked examples

```mathematica
In[1]:= GraphDiameter[CycleGraph[6]]  (* the farthest pair on a six-cycle is three steps apart *)
```

```mathematica
In[1]:= GraphDiameter[PathGraph[{1, 2, 3, 4, 5}]]  (* a path of five vertices spans four edges *)
```

```mathematica
In[1]:= GraphDiameter[GridGraph[{3, 4}]]  (* corner to corner is (3-1) + (4-1) steps *)
```

```mathematica
In[1]:= GraphDiameter[Graph[{1, 2, 3}, {UndirectedEdge[1, 2]}]]  (* an isolated vertex makes the graph disconnected *)
```

```mathematica
In[1]:= GraphDiameter[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3], DirectedEdge[3, 1]}]]  (* distances follow edge direction on a directed cycle *)
```

### Notes

The diameter is the largest shortest-path distance over all ordered pairs of vertices. An unweighted graph gives an exact integer; a graph carrying `EdgeWeight` gives a machine real.

If some vertex cannot be reached from another, the diameter is `Infinity`. For a directed graph this means the graph is not strongly connected, which is checked first in linear time so the disconnected case costs no all-pairs work.
