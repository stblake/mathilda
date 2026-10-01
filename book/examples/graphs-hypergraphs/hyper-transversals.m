# Transversals (hitting sets) and random hypergraphs.
f = Hypergraph[{{1, 2, 3}, {2, 3, 4}, {3, 4, 5}, {1, 5}}];
EdgeList[TransversalHypergraph[f]]
FindMinimumTransversal[f]
Sort[EdgeList[TransversalHypergraph[TransversalHypergraph[f]]]] == Sort[EdgeList[f]]
SeedRandom[7]; r = RandomHypergraph[{8, 5}, 3]
EdgeList[r]
EdgeList[Subhypergraph[h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}], {3, 4, 5, 6}]]
EdgeList[HypergraphRestriction[h, {3, 4, 5, 6}]]
