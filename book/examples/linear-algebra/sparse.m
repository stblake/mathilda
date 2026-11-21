# 4.4.x Sparse matrices: SparseArray, sparse Dot, direct/Krylov solvers, Krylov eigensolvers
t = SparseArray[{Band[{1, 1}] -> 2, Band[{1, 2}] -> -1, Band[{2, 1}] -> -1}, {6, 6}]
Normal[t]
Normal[t . t]
LinearSolve[t, ConstantArray[1, 6]]
Det[t]
k = 200;
t2 = SparseArray[{Band[{1, 1}] -> 2., Band[{1, 2}] -> -1., Band[{2, 1}] -> -1.}, {k, k}];
id = IdentityMatrix[k, SparseArray]
lap = KroneckerProduct[t2, id] + KroneckerProduct[id, t2]
lap["Density"]
b = ConstantArray[1., k^2];
u = LinearSolve[lap, b];
Norm[lap . u - b]
f = LinearSolve[lap]
Norm[f[b] - u]
uk = LinearSolve[lap, b, Method -> {"Krylov", Method -> "ConjugateGradient"}];
Norm[uk - u] / Norm[u]
Eigenvalues[lap, -3]
N[4 Sin[{2, 2, 1} Pi / (2 (k + 1))]^2 + 4 Sin[Pi / (2 (k + 1))]^2]
SingularValueList[lap, -2]
s = SparseArray[{Band[{1, 1}] -> 2, Band[{1, 2}] -> -1, Band[{2, 1}] -> -1}, {1000, 1000}];
LinearSolve[s, ConstantArray[1, 1000]][[500]]
LinearSolve[N[s, 40], N[ConstantArray[1, 1000], 40]][[500]]
