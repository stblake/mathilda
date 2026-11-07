(* D5 reference, method A.
   min_c ||B(c)||_2,  B(c) = A^8 + Sum_{j<8} c_j A^j,  A the 48x48 Toeplitz matrix.
   Structure used: A is Toeplitz, so J A J = A^T (J = exchange matrix), hence
   J B J = B^T and M(c) = B(c).J is SYMMETRIC and affine in c, with
   eigenvalues = +-(singular values of B).  So ||B||_2 = max(lmax(M), -lmin(M)).
   Smoothing continuation in numpy (d5_ref.py) shows that at the minimiser the
   top singular value is double (1766.31353..., next 1756.80): one copy is
   lmax(M), the other is -lmin(M) -- the two eigenvalues of opposite symmetry
   class that a generic one-parameter crossing makes equal.  Both are simple
   eigenvalues of M, hence smooth in c, and the KKT conditions of the convex
   min-max are the 9x9 smooth system
       lmax(M(c)) + lmin(M(c)) = 0,
       w grad lmax + (1 - w) grad(-lmin) = 0         (8 equations),
   grad_j l = q^T (A^j J) q (Hellmann-Feynman), unknowns (c, w), with
   0 < w < 1 as the dual-feasibility check.  Newton, Jacobian by central
   differences of the analytic gradients (h = 10^-70), fixed 250-digit
   arithmetic. *)
prec = 250;
$MinPrecision = prec; $MaxPrecision = prec;
n = 48;
A = Table[Which[i == j, 1, j == i - 1, -1, 1 <= j - i <= 3, 1, True, 0], {i, n}, {j, n}];
J = Reverse[IdentityMatrix[n]];
Mj = Table[MatrixPower[A, j].J, {j, 0, 8}];
Mof[c_] := Mj[[9]] + Sum[c[[j + 1]] Mj[[j + 1]], {j, 0, 7}];
ext[c_] := Module[{es = Eigensystem[Mof[c]], imax, imin},
   imax = First[Ordering[es[[1]], -1]]; imin = First[Ordering[es[[1]], 1]];
   {es[[1, imax]], es[[2, imax]], es[[1, imin]], es[[2, imin]]}];
F[x_] := Module[{c = x[[1 ;; 8]], w = x[[9]], lM, qM, lm, qm, gM, gm},
   {lM, qM, lm, qm} = ext[c];
   gM = Table[qM.Mj[[j]].qM, {j, 1, 8}];
   gm = Table[-qm.Mj[[j]].qm, {j, 1, 8}];
   Prepend[w gM + (1 - w) gm, lM + lm]];
c0 = SetPrecision[{1271.989218358329, -1722.6853689309014, 1464.4511145694019,
    -897.4061651604137, 419.0592230733505, -150.5652699709952,
    41.33541438702171, -7.9030639516360175}, prec];
x = Append[c0, SetPrecision[1/2, prec]];
h = SetPrecision[10^-70, prec];
Do[
  f0 = F[x];
  Print["step ", it, "  |F| = ", ScientificForm[N[Norm[f0], 20], 5],
        "  ||B|| = ", NumberForm[N[ext[x[[1 ;; 8]]][[1]], 40], 40], "  w = ", N[x[[9]], 20]];
  If[Norm[f0] < 10^-200, Break[]];
  Jac = Transpose[Table[(F[x + h UnitVector[9, k]] - F[x - h UnitVector[9, k]])/(2 h), {k, 9}]];
  x = x - LinearSolve[Jac, f0], {it, 15}];
Module[{sv = Reverse[Sort[Abs[Eigenvalues[Mof[x[[1 ;; 8]]]]]]]},
  Print["top singular values: ", N[Take[sv, 4], 30]]];
Print["c* = ", N[x[[1 ;; 8]], 30]];
Print["min ||p(A)||_2 = ", NumberForm[ext[x[[1 ;; 8]]][[1]], 80]];
