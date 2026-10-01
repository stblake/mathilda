# A Hamiltonian cycle of the 3-cube.
cube = HypercubeGraph[3];
fig = GraphPlot[cube, GraphHighlight -> First[FindHamiltonianCycle[cube]]]
