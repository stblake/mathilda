### Worked examples

```mathematica
In[1]:= HITSCentrality[Graph[{1 -> 2, 2 -> 3, 1 -> 3, 3 -> 4}]]  (* {hubs, authorities} *)
```

### Notes

The result is the pair `{h, a}`: `h` the hub scores and `a` the authority scores, each a vector
in `VertexList` order. A vertex is a good hub when it points to good authorities, and a good
authority when it is pointed to by good hubs.

`EdgeWeight` is ignored — HITS is a purely structural measure. The empty graph returns `{}`, and
a graph whose vertices never share an in-neighbour (so every co-citation class is a singleton)
scores all zeros.
