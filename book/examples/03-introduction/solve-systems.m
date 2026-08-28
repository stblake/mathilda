# Section 3.6.2: solving systems of equations -- one problem, several tools.
# A linear system: two lines meeting at a point.
LinearSolve[{{1, 1}, {1, -1}}, {3, 1}]
RowReduce[{{1, 1, 3}, {1, -1, 1}}]
# A polynomial system: the unit circle meets the line y == x.
Solve[{x^2 + y^2 == 1, y == x}, {x, y}]
GroebnerBasis[{x^2 + y^2 - 1, y - x}, {x, y}]
Resultant[x^2 + y^2 - 1, y - x, y]
# The complete solution set, as a logical description.
Reduce[{x^2 + y^2 == 1, y == x}, {x, y}]
# A system of simultaneous congruences.
ChineseRemainder[{2, 3, 2}, {3, 5, 7}]
