Solve[{x^2 + y^2 == 5, x - y == 1}, {x, y}]
Solve[{x + y + z == 6, x y + y z + z x == 11, x y z == 6}, {x, y, z}]
GroebnerBasis[{x + y + z - 6, x y + y z + z x - 11, x y z - 6}, {x, y, z}]
