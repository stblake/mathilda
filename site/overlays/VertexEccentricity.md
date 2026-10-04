### Worked examples

```mathematica
In[1]:= VertexEccentricity[PathGraph[{1, 2, 3, 4, 5}], 1]  (* an end is four hops from the far end *)
```

```mathematica
In[1]:= VertexEccentricity[PathGraph[{1, 2, 3, 4, 5}], 3]  (* the middle vertex is closer to everything *)
```

```mathematica
In[1]:= VertexEccentricity[CycleGraph[5], 1]  (* the farthest vertex on a 5-cycle is two hops away *)
```

### Notes

The eccentricity of a vertex is the greatest shortest-path distance from it to
any other vertex. The minimum eccentricity over the graph is the radius, the
maximum is the diameter, and the vertices achieving them form the center and
periphery.

For an unweighted graph the answer is an integer. With `EdgeWeight` the weights
are lengths, and a vertex that cannot reach every other has eccentricity
`Infinity`.
