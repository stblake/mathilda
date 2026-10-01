# Constraints in the global search: feasibility, infeasibility, disjunctions, far regions.
NMinimize[{x^2 + y^2, x + y >= 2}, {x, y}]
NMinimize[{x^2 + y^2, x + y == 2}, {x, y}]
NMinimize[{x + y, x^2 + y^2 <= 1 && x + y >= 3}, {x, y}]
NMinimize[{x^2, x <= -2 || x >= 2}, x]
NMinimize[{(x^2 + y - 11)^2 + (x + y^2 - 7)^2, (x + 2.8)^2 + (y + 3.1)^2 <= 0.1 || (x - 3.6)^2 + (y + 1.8)^2 <= 0.1}, {x, y}]
NMinimize[{(x - 50)^2 + (y - 40)^2, x + y >= 80}, {x, y}]
