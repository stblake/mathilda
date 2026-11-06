(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 154 -- ISSAC '97 holdout B, family 9: integral-operator and
   discretised eigenproblems (the counterpart of Problem 9, the largest
   eigenvalue of the smooth kernel E^(x+y+x^2+y^2+xy+x^2y^2) on [0, 1])
   ==========================================================================

       cd benchmarks/154-issac97b-operator-eigenvalues
       ../../Mathilda -file operator_eigenvalues.m

   The original is one smooth positive kernel on a finite interval, Nystrom
   with Gauss nodes.  These differ in kind: a kernel with a kink on the
   diagonal (algebraic Nystrom convergence, so an extrapolation or the
   secular equation), a kernel on the whole real line (Gauss-Hermite nodes),
   a periodic kernel (trapezoid nodes; a 7th eigenvalue of multiplicity 2),
   an integer matrix with a closed-form spectrum, the nonsymmetric Kac
   matrix whose integer eigenvalues are notoriously ill-conditioned, the
   Mathieu operator in a Fourier basis (and its built-in characteristic
   value), the quartic oscillator in a harmonic-oscillator basis, the sinc
   (prolate) kernel's 4th eigenvalue, the 2-D discrete Laplacian, and a
   finite-element generalized eigenproblem K x = lambda M x.
   True values (closed form, or two independent mpmath computations):

   b9.1  K = e^(-|x-y|) on [0, 1]: lambda = 2/(1 + w^2), w the least positive
         root of 2 cos w = (w - 1/w) sin w (from (1 - d^2/dx^2) K = 2 delta,
         with f' = f at 0 and f' = -f at 1, derived by hand).
         w = 1.30654237418880620222872783, lambda = 0.738810809416454997194171854902.
         Second: trapezoid Nystrom at 100/200/400 panels + two Richardson
         steps (IEEE double, numpy): 0.73881080941645, 14 digits.
   b9.2  Mehler kernel on R, rho = 1/2:
         K = exp(-((1 + rho^2)(x^2 + y^2) - 4 rho x y)/(2(1 - rho^2)))/sqrt(pi(1 - rho^2))
         = sum rho^n psi_n(x) psi_n(y) (Hermite functions; Mehler's formula),
         so the eigenvalues are 2^-n; the 6th largest is 1/32.  Second:
         Gauss-Hermite Nystrom in mpmath, n = 30 and 50 nodes (20 digits).
   b9.3  Poisson kernel (1 - r^2)/(2 pi (1 - 2 r cos(s - t) + r^2)) on the
         circle, r = 9/10: eigenvalues r^|k|, the 7th largest r^3 = 0.729.
         The N-point trapezoid discretisation has the exact eigenvalues
         (r^k + r^(N-k))/(1 - r^N) (circulant), so N = 256 is 0.729 to 4e-12
         (checked to 9 digits at machine precision, 10 at 30 digits).
   b9.4  min(i, j), N x N: lambda_max = 1/(4 sin^2(pi/(2(2N+1)))) (the inverse is
         the tridiagonal second-difference matrix), N = 500:
         101524.0106641804639020873723; N = 150: 9179.88389296667484497174920.
   b9.5  Kac (Clement) matrix, super-diagonal 1..n, sub-diagonal n..1: the
         spectrum is -n, -n+2, ..., n exactly (Kac 1947); the largest real
         part, n = 20 (machine) and n = 60 (60 digits; mpmath eig at 60
         digits returns 60.000...).
   b9.6  Mathieu a_0(q = 25) = -40.2567795465667872758017193249: the even
         pi-periodic Fourier matrix diag (2k)^2, off-diagonals q (sqrt2 q in
         the first row), truncated at 30 and 50 (identical to 30 digits);
         second: the continued-fraction characteristic equation.
   b9.7  -psi'' + x^4 psi = E psi: E_0 = 1.06036209048418289964704601669
         (published, e.g. Hioe-Montroll; mpmath diagonalisation in a
         harmonic-oscillator basis of frequency 3, sizes 60 and 100).
   b9.8  sinc kernel sin(4(x - y))/(pi (x - y)) on [-1, 1]: the 4th eigenvalue
         lambda_3(c = 4) = 0.110210987014802311973928343096 (Gauss-Legendre
         Nystrom, 30 and 50 nodes, identical to 30 digits; composite 4 x 16-point
         Gauss in IEEE double agrees to 8; lambda_0 = 0.99588549
         agrees with Slepian-Sonnenblick's table).
   b9.9  5-point Dirichlet Laplacian on an N x N grid: smallest eigenvalue
         8 sin^2(pi/(2(N+1))); N = 40: 0.0117367952650381514214343583;
         N = 10: 0.162028105542010440438527771735.
   b9.10 P1 finite elements for -u'' = lambda u on [0, pi], N interior nodes,
         h = pi/(N+1), K = tridiag(-1, 2, -1)/h, M = h tridiag(1, 4, 1)/6:
         lambda_1 = (6/h^2)(1 - cos h)/(2 + cos h); N = 400:
         1.000005114823391145999099641605; N = 60: 1.00022105340798863220108026172.

   METHODS.  Eigenvalues of the discretised operator at machine precision
   and at 30-60 digits (b9.1 also FindRoot of the secular equation; b9.6
   also MathieuCharacteristicA; b9.10 the pencil form Eigenvalues[{K, M}]).
   Node rules are computed in-file by Golub-Welsch.
   CHECKS.  Machine 6-12 digits (inside the agreement an IEEE run reaches);
   15-20 at arbitrary precision.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Eigenvalues", "Eigensystem", "FindRoot", "MathieuCharacteristicA"}];

mc[d_][v_] := If[NumericQ[v], sig[v, d], "unevaluated"];
kth[k_][e_] := Sort[Re[e], Greater][[k]];
(* symmetric Nystrom matrix sqrt(w_i) K(x_i, x_j) sqrt(w_j) *)
nys[k_, xs_, ws_] := Module[{s = Sqrt[ws]}, Table[s[[i]] k[xs[[i]], xs[[j]]] s[[j]], {i, Length[xs]}, {j, Length[xs]}]];
(* Gauss rules by Golub-Welsch: Jacobi matrix off-diagonals b(k), weight integral mu0 *)
golubWelsch[b_, n_, mu0_, prec_] := Module[{J, ev, evec},
  J = N[Table[Which[j == i + 1, b[i], i == j + 1, b[j], True, 0], {i, n}, {j, n}], prec];
  {ev, evec} = Eigensystem[J];
  {ev, mu0 evec[[All, 1]]^2}];

(* ---- b9.1 exponential kernel with a kink ------------------------------ *)
expK[x_, y_] := Exp[-Abs[x - y]];
trapNys[n_] := Module[{h = 1/n, xs, ws},
  xs = N[Range[0, n] h]; ws = N[Table[If[i == 0 || i == n, h/2, h], {i, 0, n}]];
  Max[Eigenvalues[nys[expK, xs, ws]]]];
richardson3[] := Module[{v = trapNys /@ {100, 200, 400}, r1, r2},
  r1 = (4 v[[2]] - v[[1]])/3; r2 = (4 v[[3]] - v[[2]])/3; (16 r2 - r1)/15];
run["b9.1 secular FindRoot wp40",
    2/(1 + ruleVal[FindRoot[2 Cos[w] - (w - 1/w) Sin[w], {w, 13/10}, WorkingPrecision -> 40], w]^2), mc[18], 3];
run["b9.1 Nystrom Richardson machine", richardson3[], mc[8]];

(* ---- b9.2 Mehler kernel on the line ----------------------------------- *)
mehler[x_, y_] := Exp[-(5 (x^2 + y^2) - 8 x y)/6] 2/Sqrt[3 Pi];
hermNys[n_, prec_] := Module[{xs, ws},
  {xs, ws} = golubWelsch[Sqrt[#/2] &, n, Sqrt[Pi], prec];
  Eigenvalues[nys[mehler, xs, ws Exp[xs^2]]]];
run["b9.2 Gauss-Hermite Nystrom wp30", kth[6][hermNys[40, 30]], mc[15], 3];
run["b9.2 Gauss-Hermite Nystrom machine", kth[6][hermNys[40, MachinePrecision]], mc[10], 3];

(* ---- b9.3 Poisson kernel on the circle -------------------------------- *)
poisson[s_, t_] := (1 - 81/100)/(2 Pi (1 - 18/10 Cos[s - t] + 81/100));
circNys[n_, prec_] := Module[{xs = N[2 Pi Range[0, n - 1]/n, prec]},
  Eigenvalues[nys[poisson, xs, Table[N[2 Pi/n, prec], {n}]]]];
run["b9.3 trapezoid Nystrom wp30 N=256", kth[7][circNys[256, 30]], mc[10]];
run["b9.3 trapezoid Nystrom machine N=256", kth[7][circNys[256, MachinePrecision]], mc[9], 3];

(* ---- b9.4 min(i, j) --------------------------------------------------- *)
minM[n_] := Table[Min[i, j], {i, n}, {j, n}];
run["b9.4 Eigenvalues wp50 N=150", Max[Eigenvalues[N[minM[150], 50]]], mc[20]];
run["b9.4 Eigenvalues machine N=500", Max[Eigenvalues[N[minM[500]]]], mc[10], 3];

(* ---- b9.5 Kac matrix -------------------------------------------------- *)
kac[n_] := Table[Which[j == i + 1, i, i == j + 1, n + 1 - j, True, 0], {i, n + 1}, {j, n + 1}];
run["b9.5 Eigenvalues wp60 n=60", Max[Re[Eigenvalues[N[kac[60], 60]]]], mc[15]];
run["b9.5 Eigenvalues machine n=20", Max[Re[Eigenvalues[N[kac[20]]]]], mc[6], 3];

(* ---- b9.6 Mathieu ----------------------------------------------------- *)
mathieuM[n_, q_] := Table[Which[i == j, (2 (i - 1))^2, Abs[i - j] == 1, If[Min[i, j] == 1, Sqrt[2] q, q], True, 0],
  {i, n}, {j, n}];
run["b9.6 Fourier matrix wp40", Min[Eigenvalues[N[mathieuM[40, 25], 40]]], mc[18], 3];
run["b9.6 Fourier matrix machine", Min[Eigenvalues[N[mathieuM[30, 25]]]], mc[12], 3];
run["b9.6 MathieuCharacteristicA wp30", N[MathieuCharacteristicA[0, 25], 30], mc[14]];
run["b9.6 MathieuCharacteristicA machine", MathieuCharacteristicA[0, 25.], mc[12]];

(* ---- b9.7 quartic oscillator in a harmonic basis ---------------------- *)
quartH[n_, w_, prec_] := Module[{m = n + 4, X, X2, X4},
  X = N[Table[Which[j == i + 1, Sqrt[i/(2 w)], i == j + 1, Sqrt[j/(2 w)], True, 0], {i, m}, {j, m}], prec];
  X2 = X . X; X4 = X2 . X2;
  Table[If[i == j, w (2 i - 1), 0] - w^2 X2[[i, j]] + X4[[i, j]], {i, n}, {j, n}]];
run["b9.7 harmonic basis wp50 N=100", Min[Eigenvalues[quartH[100, 3, 50]]], mc[20]];
run["b9.7 harmonic basis machine N=60", Min[Eigenvalues[quartH[60, 3, MachinePrecision]]], mc[9], 3];

(* ---- b9.8 sinc kernel ------------------------------------------------- *)
sincK[x_, y_] := If[x == y, 4/Pi, Sin[4 (x - y)]/(Pi (x - y))];
legNys[n_, prec_] := Module[{xs, ws},
  {xs, ws} = golubWelsch[#/Sqrt[4 #^2 - 1] &, n, 2, prec];
  Eigenvalues[nys[sincK, xs, ws]]];
run["b9.8 Gauss-Legendre Nystrom wp40", kth[4][legNys[30, 40]], mc[18], 3];
run["b9.8 Gauss-Legendre Nystrom machine", kth[4][legNys[30, MachinePrecision]], mc[10], 3];

(* ---- b9.9 2-D Laplacian ----------------------------------------------- *)
lap2[n_] := Module[{t = Table[Which[i == j, 2, Abs[i - j] == 1, -1, True, 0], {i, n}, {j, n}], id = IdentityMatrix[n]},
  KroneckerProduct[t, id] + KroneckerProduct[id, t]];
run["b9.9 Eigenvalues wp40 N=10", Min[Eigenvalues[N[lap2[10], 40]]], mc[18], 3];
run["b9.9 Eigenvalues machine N=40", Min[Eigenvalues[N[lap2[40]]]], mc[8]];

(* ---- b9.10 finite-element pencil -------------------------------------- *)
femKM[n_, prec_] := Module[{h = N[Pi/(n + 1), prec]},
  {Table[Which[i == j, 2/h, Abs[i - j] == 1, -1/h, True, 0], {i, n}, {j, n}],
   Table[Which[i == j, 4 h/6, Abs[i - j] == 1, h/6, True, 0], {i, n}, {j, n}]}];
run["b9.10 Eigenvalues pencil wp40 N=60", Min[Eigenvalues[femKM[60, 40]]], mc[18], 3];
run["b9.10 Eigenvalues pencil machine N=400", Min[Eigenvalues[femKM[400, MachinePrecision]]], mc[8], 3];
