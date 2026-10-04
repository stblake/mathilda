### Worked examples

```mathematica
In[1]:= GraphCenter[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5}]]  (* the one middle vertex of an odd path *)
```

### Notes

The centre is the set of vertices of minimum eccentricity — those whose greatest distance to any
other vertex equals the graph radius. It is returned as a list of vertices in `VertexList`
order.

It is `{}` unless the graph is connected (strongly connected, for a directed graph). The
matching minimum eccentricity value is `GraphRadius`.
