# Integer and mixed-integer problems: a pure integer program and its linear relaxation,
# a knapsack, and a mixed problem with an equality.
NMaximize[{x + y, -2 x + 2 y >= 1 && -8 x + 10 y <= 13 && x >= 0 && y >= 0}, {x, y}]
NMaximize[{x + y, -2 x + 2 y >= 1 && -8 x + 10 y <= 13 && x >= 0 && y >= 0 && Element[{x, y}, Integers]}, {x, y}]
NMaximize[{60 a + 100 b + 120 c, 10 a + 20 b + 30 c <= 50 && 0 <= a <= 1 && 0 <= b <= 1 && 0 <= c <= 1 && Element[{a, b, c}, Integers]}, {a, b, c}]
NMaximize[{60 a + 100 b + 120 c, 10 a + 20 b + 30 c <= 50 && 0 <= a <= 1 && 0 <= b <= 1 && 0 <= c <= 1}, {a, b, c}]
NMinimize[{x + 2 y, x^2 + 2 y^2 <= 3 && x + y == 2 && Element[x, Integers]}, {x, y}]
NMinimize[{(x - 15)^2 + (y - 3)^2, Element[y, Integers]}, {x, y}]
