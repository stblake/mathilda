### Worked examples

```mathematica
In[1]:= EdgeCount[HararyGraph[3, 8]]  (* ceil(k n / 2) = 12, the minimum for 3-connectivity *)
```

```mathematica
In[1]:= VertexConnectivity[HararyGraph[4, 7]]  (* connectivity is exactly k *)
```

```mathematica
In[1]:= VertexCount[HararyGraph[4, 7]]  (* n vertices, labelled 1..n *)
```

```mathematica
In[1]:= VertexDegree[HararyGraph[3, 8]]  (* nearly regular of degree k *)
```

### Notes

`HararyGraph[k, n]` is the extremal graph of Harary's theorem: the `k`-connected
graph on `n` vertices with as few edges as possible, namely `⌈kn/2⌉`. It is built
as the circulant on the `⌊k/2⌋` nearest neighbours of a cycle, plus one extra ring
of chords when `k` is odd. The argument order is connectivity first, so
`HararyGraph[4, 7]` is 4-connected on 7 vertices.

The result is an opaque `Graph` object; query a property (`EdgeCount`,
`VertexConnectivity`, `VertexDegree`) to inspect it. The construction requires
`n > k`, since a `k`-connected graph needs at least `k + 1` vertices; smaller `n`
is left unevaluated.
