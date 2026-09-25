# A small weighted road network among five towns; the weights are distances.
towns = {a, b, c, d, e};
links = {a <-> b, a <-> c, b <-> c, b <-> d, c <-> d, d <-> e, c <-> e};
roads = Graph[towns, links, EdgeWeight -> {2, 5, 1, 4, 2, 3, 8}];
FindShortestPath[roads, a, e]
GraphDistance[roads, a, e]
plain = Graph[towns, links];
FindShortestPath[plain, a, e]
GraphDistance[plain, a, e]
GraphDistance[Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}, EdgeWeight -> {1/2, 1/3}], 1, 3]
