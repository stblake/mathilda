### Worked examples

```mathematica
In[1]:= KatzCentrality[CycleGraph[4], 1/4]  (* a vertex-transitive graph: all scores equal *)
```

```mathematica
In[1]:= KatzCentrality[Graph[{1 -> 2, 1 -> 3}], 1/2]  (* with attenuation 1/2, each sink inherits half its source's score *)
```

```mathematica
In[1]:= KatzCentrality[Graph[{1 -> 2, 2 -> 3, 3 -> 1}], 1/10]  (* a directed triangle, small attenuation *)
```

```mathematica
In[1]:= KatzCentrality[PathGraph[{1, 2, 3}], 1/2, 2]  (* a uniform baseline b = 2 in place of 1 *)
```

### Notes

The score vector solves `x = α Aᵀ x + b`: each vertex gets a baseline `b`
(the third argument, defaulting to `1`) plus `α` times the summed scores of the
vertices pointing at it. The attenuation factor `α` must be smaller than the
reciprocal of the largest eigenvalue of the adjacency matrix for the sum to
converge; at that boundary the system is singular and the call is left
unevaluated.

Because the propagation uses `Aᵀ`, a directed edge `u -> v` lets `u`'s status flow
into `v`. Edge weights are ignored — only the graph's adjacency enters. The
baseline may also be a length-`n` list to weight the vertices differently.
