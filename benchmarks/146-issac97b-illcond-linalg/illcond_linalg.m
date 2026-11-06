(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 146 -- ISSAC '97 holdout B, family 1: ill-conditioned linear
   algebra (the counterpart of Problem 1, the Hilbert condition number)
   ==========================================================================

       cd benchmarks/146-issac97b-illcond-linalg
       ../../Mathilda -file illcond_linalg.m

   Ten problems on ten matrix families other than Hilbert, each asking a
   different ill-conditioned question (a ratio of extreme eigen- or singular
   values, the solution of a system, a least-squares fit, a determinant, a
   nearly-degenerate eigenvalue pair).  Every problem has a machine case (a
   stated reduced size where a double cannot carry the full one) and an
   arbitrary-precision case; the exact families (rational entries) also have
   an exact case.  True answers, with the two independent methods that
   established each (none of them either benchmarked system):

   b1.1  Symmetric Pascal P_n, P_ij = Binomial[i+j-2, i-1]: kappa_2.
         n = 10: 4155205697.177144144772224;  n = 60: 1.05270543880198426471473e69.
         (a) all eigenvalues at 200 digits (mpmath eigsy, Jacobi/QR);
         (b) power iteration for lambda_max alone, kappa = lambda_max^2
             (P_n is similar to its inverse, so its spectrum is closed under
             lambda -> 1/lambda).
   b1.2  Lotkin matrix (Hilbert with first row of ones), nonsymmetric:
         kappa_2 = sigma_max/sigma_min.  n = 8: 23839682200.15309347860841;
         n = 40: 2.370457593474166087334342e59.  (a) mpmath SVD at 200 digits;
         (b) square root of the eigenvalue ratio of A^T A.
   b1.3  Vandermonde system V c = b on x_i = (i-1)/(n-1), b = Runge
         1/(1 + 25 x^2): the top coefficient c_n.  n = 10:
         4010641649702490234375/26504555625524424968 = 151.3189546117181806720;
         n = 30: a 107/105-digit rational = 121163222.3328001737527216.
         (a) exact Gaussian elimination in Python Fractions; (b) mpmath LU at
         300 digits.
   b1.4  Kahan matrix R = diag(s^(i-1)) (I - c U), s = 3/5, c = 4/5, U the
         strict upper triangle of ones: sigma_min.  n = 12:
         9.956864726388756...e-6 (see below); n = 100:
         1.026723597609272540268764e-47 (sigma_max 9.905).  (a) mpmath SVD at
         160 digits; (b) the exact inverse by Fraction back-substitution, then
         power iteration on R^-T R^-1.
   b1.5  Frank matrix F_ij = n + 1 - max(i, j) for j >= i - 1 (upper
         Hessenberg, nonsymmetric; its small eigenvalues are the classic
         ill-conditioned ones): the smallest eigenvalue.  n = 12:
         0.03102806064401001507954023; n = 40: 0.007509065981719966651839415.
         (a) mpmath eig at 250 digits; (b) the exact integer characteristic
         polynomial (Hessenberg recurrence), its roots, and an exact sign
         change of the integer polynomial across the root.
   b1.6  Least-squares fit of a degree-d polynomial (monomial basis) to Runge
         at m equispaced points on [-1, 1]: the top coefficient.  m = 21,
         d = 8: 11.90607297468943863840447; m = 101, d = 20:
         4806.764449233890721337973.  (a) exact normal equations in Fractions;
         (b) mpmath Householder QR at 200 digits.
   b1.7  Hankel matrix of Catalan numbers, H_ij = C_{i+j-2}: det H = 1 for
         every n (a theorem: Hankel determinants of the Catalan sequence), and
         an integer Bareiss determinant confirms n = 7, 25, 40.  cond_2 grows
         from 3.8e9 (n = 7) to beyond 1e30.
   b1.8  Prolate matrix A_ij = sin(2 pi w (i - j))/(pi (i - j)), A_ii = 2w,
         w = 1/4 (transcendental entries, eigenvalues clustered at 0 and 1):
         the smallest eigenvalue.  n = 12: 1.778272859770121495153843e-8;
         n = 40: 1.24903564200568087964963e-29.  (a) mpmath eigsy at 150
         digits; (b) bisection on the inertia (Sylvester's law, LDL^T of
         A - sigma I).
   b1.9  Wilkinson W_{2m+1}^+ (diag |m - i|, unit off-diagonals): the gap
         between its two largest eigenvalues, a near-double pair.  m = 5:
         7.428822929311607426775169e-5; m = 20: 1.272072120490893356864979e-37
         (lambda_1 = 20.746194182903357570586883968).  (a) mpmath eigsy at 120
         digits; (b) Sturm-count bisection of the tridiagonal matrix.
   b1.10 Cauchy matrix C_ij = 1/(i^2 + j): det C.  n = 6:
         8.541575038817563368043518e-21; n = 30: a rational with a 211-digit
         numerator and an 866-digit denominator, 2.857745163831741209353925e-655.
         (a) the Cauchy determinant formula
             prod_{i<j} (x_j - x_i)(y_j - y_i) / prod_{i,j} (x_i + y_j);
         (b) exact Bareiss elimination; identical rationals.

   CHECKS.  sig[v, d] (leading digits + exponent) with d below what the method
   guarantees: machine 4-8 digits, one fewer where the true value sits within 0.05 units of a rounding edge (each machine case was also run in IEEE
   double, and d is at least two digits inside the agreement obtained), 20 at
   arbitrary precision (where at least 50 digits survive the conditioning).
   Exact rationals check {sig[v, 15], numerator mod 10^9, denominator mod 10^9}.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Eigenvalues", "SingularValueList", "LinearSolve", "LeastSquares", "Det"}];

ratChk[v_] := If[IntegerQ[v] || Head[v] === Rational,
  {sig[v, 15], Mod[Numerator[v], 10^9], Mod[Denominator[v], 10^9]}, "unevaluated"];
numQ[v_] := NumericQ[v] && FreeQ[v, Complex];

(* ---- b1.1 symmetric Pascal ------------------------------------------- *)
pascal[n_] := Table[Binomial[i + j - 2, i - 1], {i, n}, {j, n}];
kEig[m_] := Module[{e = Eigenvalues[m]}, Max[e]/Min[e]];
run["b1.1 Eigenvalues machine", kEig[N[pascal[10]]], sig[#, 5] &, 3];
run["b1.1 Eigenvalues wp200", kEig[N[pascal[60], 200]], sig[#, 20] &];

(* ---- b1.2 Lotkin ------------------------------------------------------ *)
lotkin[n_] := Table[If[i == 1, 1, 1/(i + j - 1)], {i, n}, {j, n}];
kSvd[m_] := Module[{s = SingularValueList[m]}, Max[s]/Min[s]];
run["b1.2 SingularValueList machine", kSvd[N[lotkin[8]]], sig[#, 5] &, 3];
run["b1.2 SingularValueList wp150", kSvd[N[lotkin[40], 150]], sig[#, 20] &];

(* ---- b1.3 Vandermonde solve ------------------------------------------- *)
vdm[n_] := Table[If[j == 1, 1, ((i - 1)/(n - 1))^(j - 1)], {i, n}, {j, n}];   (* 0^0 = 1 *)
runge[n_] := Table[1/(1 + 25 ((i - 1)/(n - 1))^2), {i, n}];
run["b1.3 LinearSolve exact", Last[LinearSolve[vdm[30], runge[30]]], ratChk, 3];
run["b1.3 LinearSolve wp100", Last[LinearSolve[N[vdm[30], 100], N[runge[30], 100]]], sig[#, 20] &, 3];
run["b1.3 LinearSolve machine", Last[LinearSolve[N[vdm[10]], N[runge[10]]]], sig[#, 6] &, 3];

(* ---- b1.4 Kahan ------------------------------------------------------- *)
kahan[n_] := Table[(3/5)^(i - 1) Which[j == i, 1, j > i, -4/5, True, 0], {i, n}, {j, n}];
run["b1.4 SingularValueList machine", Min[SingularValueList[N[kahan[12]]]], sig[#, 5] &, 3];
run["b1.4 SingularValueList wp100", Min[SingularValueList[N[kahan[100], 100]]], sig[#, 20] &];

(* ---- b1.5 Frank ------------------------------------------------------- *)
frank[n_] := Table[If[j >= i - 1, n + 1 - Max[i, j], 0], {i, n}, {j, n}];
minAbsEig[m_] := Module[{e = Eigenvalues[m]}, Re[First[SortBy[e, Abs]]]];
run["b1.5 Eigenvalues machine", minAbsEig[N[frank[12]]], sig[#, 4] &, 3];
run["b1.5 Eigenvalues wp150", minAbsEig[N[frank[40], 150]], sig[#, 20] &];

(* ---- b1.6 polynomial least squares ------------------------------------ *)
lsA[m_, d_] := Table[If[j == 0, 1, (-1 + 2 (i - 1)/(m - 1))^j], {i, m}, {j, 0, d}];
lsB[m_] := Table[1/(1 + 25 (-1 + 2 (i - 1)/(m - 1))^2), {i, m}];
run["b1.6 LeastSquares exact", Last[LeastSquares[lsA[101, 20], lsB[101]]], ratChk, 3];
run["b1.6 LeastSquares wp60", Last[LeastSquares[N[lsA[101, 20], 60], N[lsB[101], 60]]], sig[#, 20] &, 3];
run["b1.6 LeastSquares machine", Last[LeastSquares[N[lsA[21, 8]], N[lsB[21]]]], sig[#, 8] &, 3];

(* ---- b1.7 Catalan Hankel determinant ---------------------------------- *)
catHankel[n_] := Table[Binomial[2 (i + j - 2), i + j - 2]/(i + j - 1), {i, n}, {j, n}];
oneChk[tol_][v_] := If[NumericQ[v], Abs[v - 1] < tol, "unevaluated"];
run["b1.7 Det exact", Det[catHankel[40]], If[IntegerQ[#], #, "unevaluated"] &, 3];
run["b1.7 Det wp120", Det[N[catHankel[25], 120]], oneChk[10^-20], 3];
run["b1.7 Det machine", Det[N[catHankel[7]]], oneChk[10^-6], 3];

(* ---- b1.8 prolate matrix ---------------------------------------------- *)
prolate[n_] := Table[If[i == j, 1/2, Sin[Pi (i - j)/2]/(Pi (i - j))], {i, n}, {j, n}];
run["b1.8 Eigenvalues machine", Min[Eigenvalues[N[prolate[12]]]], sig[#, 4] &, 3];
run["b1.8 Eigenvalues wp120", Min[Eigenvalues[N[prolate[40], 120]]], sig[#, 20] &];

(* ---- b1.9 Wilkinson W+ near-double pair ------------------------------- *)
wilk[m_] := Table[Which[i == j, Abs[m + 1 - i], Abs[i - j] == 1, 1, True, 0], {i, 2 m + 1}, {j, 2 m + 1}];
gap[mat_] := Module[{e = Sort[Eigenvalues[mat], Greater]}, e[[1]] - e[[2]]];
run["b1.9 Eigenvalues machine", gap[N[wilk[5]]], sig[#, 7] &, 3];
run["b1.9 Eigenvalues wp60", gap[N[wilk[20], 60]], sig[#, 15] &];

(* ---- b1.10 Cauchy determinant ----------------------------------------- *)
cauchy[n_] := Table[1/(i^2 + j), {i, n}, {j, n}];
run["b1.10 Det exact", Det[cauchy[30]], ratChk, 3];
run["b1.10 Det wp300", Det[N[cauchy[30], 300]], sig[#, 20] &, 3];
run["b1.10 Det machine", Det[N[cauchy[6]]], sig[#, 5] &, 3];
