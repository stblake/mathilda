# Flows and cuts.
net = Graph[{s, a, b, t}, {s -> a, s -> b, a -> b, a -> t, b -> t}]
FindMaximumFlow[net, s, t, EdgeCapacity -> {3, 2, 1, 2, 3}]
FindMaximumFlow[net, s, t, "FlowMatrix", EdgeCapacity -> {3, 2, 1, 2, 3}]
FindEdgeCut[net, s, t]
g = PetersenGraph[];
FindMinimumCut[g]
FindVertexCut[g]
{VertexConnectivity[g], EdgeConnectivity[g], Min[VertexDegree[g]]}
