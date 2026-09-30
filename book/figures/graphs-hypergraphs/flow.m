# A maximum flow, with each edge labelled flow/capacity.
net = Graph[{s, a, b, t}, {s -> a, s -> b, a -> b, a -> t, b -> t}];
fig = GraphPlot[net, VertexLabels -> "Name", EdgeLabels -> {(s -> a) -> "3/3", (s -> b) -> "2/2", (a -> b) -> "1/1", (a -> t) -> "2/2", (b -> t) -> "3/3"}, GraphLayout -> "LayeredDigraphEmbedding"]
