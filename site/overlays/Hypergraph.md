### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]  (* the running example: 7 vertices, 4 hyperedges *)
```

```mathematica
In[1]:= VertexList[h]  (* vertices, in first-appearance order *)
```

```mathematica
In[1]:= EdgeList[h]  (* the hyperedges, kept exactly as written *)
```

```mathematica
In[1]:= {VertexCount[h], EdgeCount[h]}  (* the shared Graph accessors work on a Hypergraph *)
```

### Notes

A `Hypergraph` is an ordinary `Expr` tree, `Hypergraph[{v1, ..., vn}, {e1, ..., em}]`,
with no new kernel type — the same uniformity that lets `Part`, `Map` and the
Graph accessors apply to it. Each hyperedge is a `List` of vertices; vertices are
arbitrary, pairwise-distinct expressions. Hyperedges may repeat (a
multi-hypergraph), overlap, nest, be empty, or repeat a vertex (`{1, 1, 2}`, as
Wolfram-model states do).

Order is kept, and that is deliberate. Order-sensitive heads
(`HypergraphToGraph`, `EdgeList`, the arity heads, `HypergraphEdgeDelete`) see the
hyperedge exactly as written; every set-theoretic head reads it as the set of its
distinct vertices. So one object answers both the ordered (directed) and the
unordered question.

Standard output is the terse summary `Hypergraph[<n vertices, m hyperedges>]`;
`InputForm` and `FullForm` print the literal, which round-trips through the parser.
Mathematica has no built-in `Hypergraph`; the name follows the
WolframInstitute/Hypergraph paclet.
