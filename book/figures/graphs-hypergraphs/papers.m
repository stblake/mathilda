# Five papers as hyperedges over their authors.
papers = Hypergraph[{{"Erdos", "Renyi"}, {"Erdos", "Graham", "Chung"}, {"Graham", "Knuth", "Patashnik"}, {"Knuth", "Yao"}, {"Chung", "Yao", "Graham"}}];
fig = HypergraphPlot[papers, VertexLabels -> "Name"]
