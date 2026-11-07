(* Holdout A, family hd5, method A (wolframscript), independent code from
   hd5_ref.py:
   - hd5.3 (Frobenius): EXACT rational normal equations, value Sqrt[rational].
   - normal matrices (hd5.1, .2, .4, .7): the discrete minimax as an exact /
     100-digit LinearOptimization (min t s.t. -t <= p(mu_i) <= t).
   - non-normal Toeplitz (hd5.5, .6, .8, .9, .10): as d5_ref.wl -- M(c) = B(c) J
     symmetric, KKT Newton on the two-fold crease lmax(M) + lmin(M) = 0,
     w grad lmax + (1-w) grad(-lmin) = 0, Eigensystem at fixed 120 digits,
     central-difference Jacobian; started from the stage-1 coefficients that
     hd5_ref.py prints (hd5_start.wl).  Run: wolframscript -file hd5_ref.wl *)
toep[n_, d_] := Table[Lookup[d, j - i, 0], {i, n}, {j, n}];
grcar[n_, s_] := toep[n, Join[<|-1 -> -1, 0 -> 1|>, AssociationThread[Range[s] -> 1]]];
basis[A_, k_, nz_] := Module[{pw = Table[MatrixPower[A, j], {j, 0, k}]},
   If[nz === "monic", {pw[[k + 1]], pw[[1 ;; k]]}, {pw[[1]], pw[[2 ;; k + 1]]}]];
(* Frobenius *)
Module[{b0, bs, G, g, c},
  {b0, bs} = basis[grcar[48, 3], 8, "monic"];
  G = Table[Total[Flatten[a b]], {a, bs}, {b, bs}]; g = Table[Total[Flatten[a b0]], {a, bs}];
  c = -LinearSolve[G, g];
  Print["3 ", NumberForm[N[Sqrt[Total[Flatten[b0^2]] + c.g], 80], 75]]];
(* normal: discrete minimax by LP on the spectrum *)
lp[mus_, k_, nz_, prec_] := Module[{pw, fx, cs, t, cons, sol, pw0},
  {pw, fx} = If[nz === "monic", {Range[0, k - 1], k}, {Range[k], 0}];
  cs = Array[cc, Length[pw]];
  pw0[m_, j_] := If[j == 0, 1, m^j];   (* 0^0 = 1: Cos[Pi/2] = 0 is an eigenvalue in hd5.1 *)
  cons = Flatten[Table[{cs.(pw0[m, #] & /@ pw) + pw0[m, fx] <= t, -t <= cs.(pw0[m, #] & /@ pw) + pw0[m, fx]}, {m, mus}]];
  sol = If[prec === Infinity,
     LinearOptimization[t, cons, Append[cs, t]],
     LinearOptimization[t, cons, Append[cs, t], WorkingPrecision -> prec]];
  t /. sol];
Print["1 ", NumberForm[lp[N[Table[Cos[j Pi/8], {j, 0, 8}], 100], 8, "monic", 100], 75]];
Print["2 ", lp[Range[10], 4, "monic", Infinity]];
Print["4 ", NumberForm[lp[N[Table[2 - 2 Cos[j Pi/21], {j, 20}], 100], 5, "monic", 100], 75]];
Print["7 ", lp[Range[10], 3, "gmres", Infinity]];
(* hd5.8: A = I + 2 N (N the nilpotent shift, n = 12 > k = 4): the Chebyshev
   polynomial of a Jordan block is (z - 1)^k (Greenbaum-Trefethen 1994), so the
   minimum is ||(2N)^4||_2 = 16 exactly; check it and that the smoothing stage of
   hd5_ref.py reaches 16.00000000000000 with an 8-fold top singular value. *)
Print["8 ", Norm[MatrixPower[toep[12, <|0 -> 1, 1 -> 2|>] - IdentityMatrix[12], 4]]];
(* non-normal: two-fold crease Newton *)
Get["hd5_start.wl"];
prec = 120; $MinPrecision = prec; $MaxPrecision = prec;
crease[A_, k_, c0_] := Module[{n = Length[A], J, b0, bs, Mj, Mof, ext, F, x, h, f0, Jac},
  J = Reverse[IdentityMatrix[n]];
  {b0, bs} = basis[A, k, "monic"];
  Mj = N[Prepend[(# . J) & /@ bs, b0 . J], prec];
  Mof[c_] := Mj[[1]] + c . Mj[[2 ;;]];
  ext[c_] := Module[{es = Eigensystem[Mof[c]], imax, imin},
    imax = First[Ordering[es[[1]], -1]]; imin = First[Ordering[es[[1]], 1]];
    {es[[1, imax]], es[[2, imax]], es[[1, imin]], es[[2, imin]]}];
  F[xx_] := Module[{c = xx[[1 ;; k]], w = xx[[k + 1]], lM, qM, lm, qm},
    {lM, qM, lm, qm} = ext[c];
    Prepend[w Table[qM . Mj[[j + 1]] . qM, {j, k}] - (1 - w) Table[qm . Mj[[j + 1]] . qm, {j, k}],
            lM + lm]];
  x = Append[SetPrecision[c0, prec], SetPrecision[1/2, prec]];
  h = SetPrecision[10^-45, prec];
  Do[f0 = F[x]; If[Norm[f0] < 10^-95, Break[]];
     Jac = Transpose[Table[(F[x + h UnitVector[k + 1, q]] - F[x - h UnitVector[k + 1, q]])/(2 h), {q, k + 1}]];
     x = x - LinearSolve[Jac, f0], {it, 30}];
  Print["   |F| = ", N[Norm[F[x]], 3], "  w = ", N[x[[k + 1]], 10],
        "  top |eig| = ", N[Take[Reverse[Sort[Abs[Eigenvalues[Mof[x[[1 ;; k]]]]]]], 3], 12]];
  ext[x[[1 ;; k]]][[1]]];
Do[With[{A = pr[[1]], k = pr[[2]], id = pr[[3]]},
   Print[id, " ", NumberForm[crease[A, k, hd5start[id]], 75]]],
  {pr, {{grcar[32, 3], 6, 5}, {grcar[40, 2], 6, 6},{toep[24, <|-1 -> 1/4, 1 -> 1|>], 6, 9}, {grcar[64, 3], 10, 10}}}];
