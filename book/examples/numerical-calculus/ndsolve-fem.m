# Finite elements: eigenvalues, Poisson's equation, and integrals over a mesh
NDEigenvalues[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, u[x, y], {x, y} \[Element] Rectangle[{0, 0}, {Pi, Pi}], 4]
NDEigenvalues[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, u[x, y], {x, y} \[Element] Disk[], 3]
N[BesselJZero[0, 1]^2]
NDEigenvalues[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, u[x, y], {x, y} \[Element] Polygon[{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}], 1]
v = NDSolveValue[{-Laplacian[u[x, y], {x, y}] == 1, DirichletCondition[u[x, y] == 0, True]}, u, {x, y} \[Element] Disk[]]
{v[0, 0], v[0.5, 0]}
NIntegrate[v[x, y], {x, y} \[Element] Disk[]]
N[Pi/8]
