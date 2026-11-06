(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 144 -- ISSAC '97 holdout, family 9: eigenvalues of operators
   ==========================================================================

       cd benchmarks/144-issac97h-operator-eigenvalues
       ../../Mathilda -file operator_eigenvalues.m

   HOLDOUT.  Ten problems probing what experiment 134 (the largest eigenvalue
   of the integral operator with kernel E^(x+y+x^2+y^2+xy+x^2y^2)) probes --
   a dominant (or second) eigenvalue of an integral operator by Nystrom
   discretisation, and eigenvalues of discretised differential operators --
   on different kernels and operators.  See benchmarks/ISSAC97_HOLDOUT.md.

   h9.1-h9.6, h9.10: Nystrom with composite 32-point Gauss-Legendre (nodes by
   Golub-Welsch, as in experiment 134) on c = 2^i panels, i increased until
   the eigenvalue is stationary to tol (machine 10^-12, checked to 9 digits;
   WorkingPrecision 40: 10^-22, checked to 18 digits).  The matrix is
   symmetrised, W^(1/2) K W^(1/2), so Eigenvalues sees a symmetric matrix.
   h9.7-h9.9: a fixed discretisation (Fourier, Hermite-function and sine
   Galerkin bases), machine and 40 digits, sized so the discretisation error
   is below the check.

   PROBLEMS AND TRUE VALUES (the Nystrom values at 60 and 90 Gauss nodes and
   80 digits agree to >= 33 digits; each second method is named):
   h9.1  lambda_1, K = E^(x y) on [0, 1]          = 1.353030164745734853309201750
   h9.2  lambda_1, K = 1/(1 + x + y) on [0, 1]    = 0.5362060346834763502347030920
   h9.3  lambda_1, K = Cos[x y] on [0, Pi]        = 1.253314101714016907331673330
   h9.4  lambda_2, K = E^(-(x - y)^2) on [0, 2]   = 0.5359573175891977563406853434
   h9.5  lambda_1, K = 1/(1 + (x - y)^2) on [-1, 1] = 1.424941177795263197428278867
   h9.6  lambda_1, K = BesselJ[0, x y] on [0, 4]  = 1.665886980965023224079412840
   h9.7  Mathieu: the smallest eigenvalue of -y'' + 2 q Cos[2 x] y on even
         Pi-periodic functions, q = 5, = MathieuCharacteristicA[0, 5]
         = -5.800046020851508419103472514 (the special function, and the
         Fourier-cosine matrix at 60 modes; machine 20 modes, wp40 40 modes)
   h9.8  ground state of -y'' + x^4 y = E y on the line
         = 1.060362090484182899647046017 (published: Hioe and Montroll;
         Rayleigh-Ritz in 100 and 200 Hermite functions of frequency 3 agree
         to 41 digits; machine 40 functions, wp40 100)
   h9.9  smallest eigenvalue of -y'' = lambda (1 + x) y, y(0) = y(1) = 0
         = 6.548395306000592515133370778 (the Airy-function characteristic
         equation Ai(-l^(1/3)) Bi(-2 l^(1/3)) = Ai(-2 l^(1/3)) Bi(-l^(1/3)),
         solved at 60 digits; the generalized sine-Galerkin eigenproblem
         A c = lambda B c converges to it algebraically -- 1e-11 at 20 modes
         (machine, 7 digits checked), 1e-17 at 200 modes (wp40, 14 digits))
   h9.10 lambda_1, K = x y + E^(x + y) on [0, 1] (rank 2)
         = (-1 + 3 E^2 + Sqrt[169 - 30 E^2 + 9 E^4])/12
         = 3.509384443552449794992539729 (closed form: the eigenvalues of the
         2 x 2 Gram matrix of {x, E^x})
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Eigensystem", "Eigenvalues"}];

(* n-point Gauss-Legendre on (0, 1) at precision prec, by Golub-Welsch. *)
gaussLegendre01[n_, prec_] := Module[{J, ev, evec, b},
  b[k_] := k/Sqrt[4 k^2 - 1];
  J = N[Table[Which[j == i + 1, b[i], i == j + 1, b[j], True, 0], {i, n}, {j, n}], prec];
  {ev, evec} = Eigensystem[J];
  Transpose[Sort[Transpose[{(1 + ev)/2, evec[[All, 1]]^2}]]]];

(* idx-th largest eigenvalue of the kernel k on [a, b]; panels doubled until
   the relative change is below tol (at most 2^4 panels, as in exp. 134). *)
nystrom[k_, a_, b_, idx_, prec_, tol_] := Module[{xs, ws, eOld = 1, eNew, i = 0, c, pts, wts, s},
  {xs, ws} = gaussLegendre01[32, prec];
  While[c = 2^i;
    pts = a + (b - a) Flatten[Table[xs + j, {j, 0, c - 1}]/c];
    wts = (b - a) Flatten[Table[ws, {c}]/c];
    s = Sqrt[wts];
    eNew = Sort[Re[Eigenvalues[Outer[Times, s, s] Outer[k, pts, pts]]], Greater][[idx]];
    Abs[(eNew - eOld)/eNew] > tol && i < 4,
    i = i + 1; eOld = eNew];
  eNew];

kernels = {
  {"h9.1", Exp[#1 #2] &, 0, 1, 1},
  {"h9.2", 1/(1 + #1 + #2) &, 0, 1, 1},
  {"h9.3", Cos[#1 #2] &, 0, Pi, 1},
  {"h9.4", Exp[-(#1 - #2)^2] &, 0, 2, 2},
  {"h9.5", 1/(1 + (#1 - #2)^2) &, -1, 1, 1},
  {"h9.6", BesselJ[0, #1 #2] &, 0, 4, 1}};

Do[With[{lab = kk[[1]], k = kk[[2]], a = kk[[3]], b = kk[[4]], i = kk[[5]]},
   run[lab <> " Nystrom machine", nystrom[k, a, b, i, MachinePrecision, 10^-12], sig[#, 9] &, 3];
   run[lab <> " Nystrom wp40", nystrom[k, a, b, i, 40, 10^-22], sig[#, 18] &]],
  {kk, kernels}];

(* h9.7 Mathieu, even Pi-periodic: basis 1/Sqrt[2], Cos[2x], Cos[4x], ... *)
mathieu[n_, q_] := Table[Which[i == j, 4 (i - 1)^2, Abs[i - j] == 1, If[Min[i, j] == 1, Sqrt[2] q, q], True, 0],
  {i, n}, {j, n}];
run["h9.7 Fourier matrix machine n=20", Min[Eigenvalues[N[mathieu[20, 5]]]], sig[#, 9] &, 3];
run["h9.7 Fourier matrix wp40 n=40", Min[Eigenvalues[N[mathieu[40, 5], 40]]], sig[#, 18] &];

(* h9.8 -d^2/dx^2 + x^4 in n Hermite functions of frequency w:
   x = (a + a^+)/Sqrt[2 w], p^2 = w (2 i - 1)/2 on the diagonal, -w Sqrt[(i-1)(i-2)]/2 at offset 2. *)
xmat[n_, w_] := Table[Which[j == i + 1, Sqrt[i]/Sqrt[2 w], i == j + 1, Sqrt[j]/Sqrt[2 w], True, 0], {i, n}, {j, n}];
anharm[n_, w_, p_] := Module[{X = N[xmat[n + 4, w], p], P2},
  P2 = N[Table[Which[i == j, w (2 i - 1)/2, Abs[i - j] == 2, -w Sqrt[Max[i, j] - 1] Sqrt[Max[i, j] - 2]/2, True, 0],
     {i, n}, {j, n}], p];
  Min[Eigenvalues[P2 + (X . X . X . X)[[1 ;; n, 1 ;; n]]]]];
run["h9.8 Hermite basis machine n=40", anharm[40, 3, MachinePrecision], sig[#, 9] &, 3];
run["h9.8 Hermite basis wp40 n=100", anharm[100, 3, 40], sig[#, 18] &];

(* h9.9 sine Galerkin: A = diag((k Pi)^2/2), B_jk = Int_0^1 (1 + x) Sin[j Pi x] Sin[k Pi x]. *)
cc[k_] := If[k == 0, 1/2, ((-1)^k - 1)/(k^2 Pi^2)];
galerkin[n_, p_] := Min[Re[Eigenvalues[{N[DiagonalMatrix[Table[(i Pi)^2/2, {i, n}]], p],
     N[Table[If[i == j, 1/2, 0] + (cc[i - j] - cc[i + j])/2, {i, n}, {j, n}], p]}]]];
run["h9.9 sine Galerkin machine n=20", galerkin[20, MachinePrecision], sig[#, 7] &, 3];
run["h9.9 sine Galerkin wp40 n=200", galerkin[200, 40], sig[#, 14] &];

run["h9.10 Nystrom machine", nystrom[#1 #2 + Exp[#1 + #2] &, 0, 1, 1, MachinePrecision, 10^-12], sig[#, 9] &, 3];
run["h9.10 Nystrom wp40", nystrom[#1 #2 + Exp[#1 + #2] &, 0, 1, 1, 40, 10^-22], sig[#, 18] &];
