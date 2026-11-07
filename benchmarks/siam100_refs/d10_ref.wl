(* Reference for d10 (random recurrence x_{n+1} = 2 x_n +- x_{n-1}), method A.
   r_n = x_{n+1}/x_n stays in [1,3] and follows the random Moebius maps
   r -> 2 - 1/r, r -> 2 + 1/r (prob 1/2 each).  The a.s. limit of |x_n|^(1/n)
   is exp(lambda), lambda = Int log r dnu(r), nu the stationary law of r
   (Furstenberg).  lambda is the constant in the Poisson equation
       h(r) - (h(2-1/r) + h(2+1/r))/2 = log r - lambda   on [1,3],
   solved by Chebyshev collocation (N nodes, h = Sum_{k>=1} c_k T_k(r-2)),
   at WorkingPrecision wp, for increasing N. *)
lam[nn_, wp_] := Module[{xs, row, A, b, sol},
  xs = Table[N[Cos[Pi (j - 1/2)/nn], wp], {j, nn}];
  row[z_] := Table[ChebyshevT[k, z], {k, 1, nn - 1}];
  A = Table[With[{r = 2 + x}, Append[row[x] - (row[-1/r] + row[1/r])/2, 1]], {x, xs}];
  b = Log[2 + xs];
  sol = LinearSolve[A, b];
  Exp[Last[sol]]];
Do[Print[nn, "  ", N[lam[nn, 120], 75]], {nn, {200, 300, 400, 500, 600, 700}}];
