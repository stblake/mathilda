### Worked examples

```mathematica
In[1]:= VertexCount[HypercubeGraph[0]]  (* the zero-dimensional cube is a single vertex *)
```

```mathematica
In[1]:= EdgeList[HypercubeGraph[1]]  (* the cube of dimension one is one edge *)
```

```mathematica
In[1]:= EdgeList[HypercubeGraph[2]]  (* the square, a 4-cycle *)
```

```mathematica
In[1]:= VertexDegree[HypercubeGraph[3]]  (* every vertex of the d-cube has degree d *)
```

```mathematica
In[1]:= EdgeCount[HypercubeGraph[4]]  (* d times 2^(d-1) edges *)
```

```mathematica
In[1]:= AdjacencyMatrix[HypercubeGraph[2]]  (* vertices adjacent when their labels differ in one bit *)
```

### Notes

`HypercubeGraph[d]` has `2^d` vertices; vertices `i` and `j` are adjacent exactly when `i - 1` and `j - 1` differ in a single binary digit. It is `d`-regular and bipartite.

The dimension must be an integer from 0 to 24; anything else leaves the call unevaluated.
