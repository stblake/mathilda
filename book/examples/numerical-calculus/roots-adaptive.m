# High-degree roots, complete solution sets, FindRoot at WorkingPrecision
r = x /. NSolve[LaguerreL[200, x] == 0, x, WorkingPrecision -> 30];
{Length[r], Count[r, _Real], Max[r]}
Length[NSolve[{x^2 + y^2 + z^2 == 1, x + y + z == 0, x y z == 1/10}, {x, y, z}]]
FindRoot[LaguerreL[1000, x], {x, 3990, 4002}, WorkingPrecision -> 40]
FindRoot[LaguerreL[1000, x], {x, 4002}]
FindRoot[Exp[1000 x] - 2, {x, 1}]
