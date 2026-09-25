# 4.4.4 Solving linear systems: LinearSolve and LeastSquares
LinearSolve[{{1, 2}, {3, 4}}, {5, 6}]
LinearSolve[{{r, s}, {t, u}}, {y, z}]
LinearSolve[{{1, 2, 3}, {4, 5, 6}}, {6, 15}]
LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}]
LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 7}] == LinearSolve[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 7}]
