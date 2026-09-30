# Many variables: indexed variables in NMinimize, generated symbols in FindMinimum.
NMinimize[Sum[x[i]^2 - 10 Cos[2 Pi x[i]] + 10, {i, 1, 5}], Table[x[i], {i, 1, 5}]]
NMinimize[Sum[100 (x[i + 1] - x[i]^2)^2 + (1 - x[i])^2, {i, 1, 4}], Array[x, 5]]
vars = Table[Symbol["z" <> ToString[i]], {i, 6}]
f6 = Sum[100 (vars[[i + 1]] - vars[[i]]^2)^2 + (1 - vars[[i]])^2, {i, 1, 5}];
FindMinimum[Evaluate[f6], Evaluate[Transpose[{vars, ConstantArray[0, 6]}]], Method -> "LBFGSB"]
