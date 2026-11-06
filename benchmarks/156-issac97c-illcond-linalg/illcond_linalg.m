(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 156 -- ISSAC '97 holdout C (SEALED), family c1:
   ill-conditioned linear algebra (counterpart of Problem 1)
   ==========================================================================

       cd benchmarks/156-issac97c-illcond-linalg
       ../../Mathilda -file illcond_linalg.m

   Ten matrix structures, none of them Hilbert-like: a Vandermonde on the
   integers, a Gaussian (theta) Toeplitz, the Moler matrix, a grounded graph
   Laplacian with geometric conductances, a near-singular non-symmetric
   circulant, Kac-Murdock-Szego with rho -> 1, Hankel matrices of the Bernoulli
   numbers and of the factorials, the Lauchli least-squares problem and the
   (non-symmetric, diagonally dominant) Dorr matrix.  Each problem has a
   machine case at the largest size a double still resolves, an arbitrary-
   precision case (the matrix is entered with N[m, wp]) and, where the answer is
   rational, an exact case.  kappa = sigma_max/sigma_min (SingularValueList) or
   lambda_max/lambda_min (Eigenvalues, symmetric positive definite matrices).

   TRUE ANSWERS (A = high-precision mpmath svd_r/eigsy/lu_solve at >= 2x the
   digits lost, or exact Fraction arithmetic; B = an independent route: a
   closed form derived by hand, or PARI's Jacobi eigenvalues qfjacobi of the
   exact rational matrix, at 80-600 digits).  The three largest wp cases
   (Vandermonde n=60, Moler n=150, factorial Hankel n=60) are beyond mpmath's
   dense SVD; there A = qfjacobi of the matrix (or of V^T V) at 700 digits and
   B = qfjacobi of its exact inverse (for Moler the closed-form inverse
   L^-T L^-1, checked M M^-1 = I), whose extreme eigenvalue ratio is the same
   kappa computed from the other end of the spectrum.  A and B agree to every
   digit shown.  Digits checked: machine 5-8 (at or below -log10(eps kappa)), wp 20
   (at most half the working precision left after the digits kappa costs).

   c1.1  V_ij = i^(j-1).  exact Det n=30 = prod_{k=1}^{29} k! (superfactorial,
         B) = 2.15551190904274323928924e383, Mod 1e9+7 = 688963510;
         machine kappa n=8 = 952111739.17589770431;
         wp300 kappa n=60 = 6.50027973410168833639888752219e123.
   c1.2  G_ij = q^((i-j)^2).  B: G = D V D with V the Vandermonde in q^(-2i),
         so Det = prod_{k=1}^{n-1} (1 - q^(2k))^(n-k).  exact Det n=40, q=9/10
         = 4.5779793077473346017e-95 (numerator Mod 1e9+7 = 680396165,
         denominator 786603587); machine lambda_min n=30, q=9/10
         = 4.6374958257459897152e-9; wp60 lambda_min n=120, q=19/20
         = 3.67554667929975235737236441273e-20.
   c1.3  Moler M = L L^T, L unit lower triangular with -1 below the diagonal
         (M_ii = i, M_ij = min(i,j) - 2).  B: L^-1 has entries 2^(i-j-1) below
         the diagonal, so (M^-1)_11 = 1 + (4^(n-1) - 1)/3.  exact
         Inverse[M_30][[1,1]] = 96076792050570582; machine kappa n=12
         = 87522898.09438388228; wp250 kappa n=150
         = 2.02335096218540765443818115297e93.
   c1.4  Grounded path Laplacian, conductance 2^-k on edge (k-1,k), node 0
         grounded.  B: a unit current into node n gives the series-resistance
         voltages x_j = sum_{k<=j} 2^k = 2^(j+1) - 2.  exact LinearSolve n=60:
         {x_60, x_30} = {2305843009213693950, 2147483646}; machine kappa n=20
         = 2682440.2336704333605; wp60 kappa n=60
         = 2949381270885395289.11123505789.
   c1.5  Circulant (1+eps) I - P (P the cyclic shift): normal, singular values
         |1 + eps - w^k|, so for even n kappa = (2 + eps)/eps (B).  machine
         n=32, eps=1e-6: 2000001; wp60 n=48, eps=1e-30: 2e30 + 1.
   c1.6  KMS A_ij = r^|i-j|.  B: Det = (1 - r^2)^(n-1) (the inverse is
         tridiagonal).  exact Det n=50, r=999/1000 = 5.4932189779635338119e-133
         (numerator Mod 1e9+7 = 892345966, denominator 335996330); machine
         lambda_min n=50, r=0.999 = 5.0074416716795090139e-4; wp60
         lambda_min n=50, r = 1 - 1e-12 = 5.00493805098963979720392692316e-13.
   c1.7  Hankel of Bernoulli numbers H_ij = B_(i+j-2) (B_1 = -1/2).  B: the
         Al-Salam-Carlitz evaluation Det = (-1)^(n(n-1)/2) prod_{k=1}^{n-1}
         k!^6/((2k)!(2k+1)!), checked equal to the Fraction determinant.  exact
         Det n=14 = -3.1050410217561563263 (numerator Mod 1e9+7 = 351202457,
         denominator 287719374); machine kappa n=8 = 7885.9266359339996083;
         wp60 kappa n=20 = 224185454513517586077.459984416.
   c1.8  Lauchli: A = [1^T; eps I_n], b = (1, 1, 2, ..., n), n = 20.  B: by
         Sherman-Morrison on A^T A = 1 1^T + eps^2 I,
         x_i = (eps + i (n + eps^2) - n(n+1)/2)/(eps (n + eps^2)).  x_20:
         exact eps=10^-10: 190000000000100000000200000000000/
         2000000000000000000001 (Mods 934833507, 98001); machine eps=1e-6:
         9500000.050000525; wp80 eps=1e-30: 9500000000000000000000000000000.05.
         The normal equations lose everything at machine precision; QR does not.
   c1.9  Hankel of factorials H_ij = (i+j-2)! (moments of e^-x).  B: Det =
         prod_{k=0}^{n-1} k!^2 (Laguerre norms).  exact Det n=30
         = 4.6462315900250914036e766 (Mod 1e9+7 = 788825102); machine kappa n=6
         = 350717286.23187149429; wp400 kappa n=60
         = 3.14466193056846718045811323198e205.
   c1.10 Dorr matrix D(n, theta) (definition below, as in the classical test
         gallery), solve D x = 1, report x_(n/2).  A: Fraction solve;
         B: mpmath lu_solve.  exact n=50, theta=1/100: 5868.37134401954916514
         (kappa 1.3e6); machine n=60, theta=1/200: 8081580.43332516446719
         (kappa 1.9e9); wp60 n=100, theta=1/1000:
         2261277808846905453744539.11567 (kappa 8.6e26).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"SingularValueList", "Eigenvalues", "LinearSolve", "LeastSquares", "Det", "Inverse", "BernoulliB"}];

exactChk[v_] := If[IntegerQ[v], {sig[v, 15], Mod[v, 1000000007]}, "unevaluated"];
ratChk[v_] := If[IntegerQ[v] || Head[v] === Rational,
  {sig[v, 15], Mod[Numerator[v], 1000000007], Mod[Denominator[v], 1000000007]}, "unevaluated"];
kappa[m_] := Module[{s = SingularValueList[m]}, Max[s]/Min[s]];
symKappa[m_] := Module[{e = Eigenvalues[m]}, Max[e]/Min[e]];
lmin[m_] := Min[Eigenvalues[m]];

vand[n_] := Table[i^(j - 1), {i, n}, {j, n}];
gtoep[n_, q_] := Table[q^((i - j)^2), {i, n}, {j, n}];
moler[n_] := Table[If[i == j, i, Min[i, j] - 2], {i, n}, {j, n}];
lap[n_] := Module[{c = Append[Table[2^-k, {k, n}], 0]},
  Table[Which[i == j, c[[i]] + c[[i + 1]], j == i + 1, -c[[i + 1]], i == j + 1, -c[[j + 1]], True, 0],
        {i, n}, {j, n}]];
circ[n_, eps_] := Table[Which[i == j, 1 + eps, j == Mod[i, n] + 1, -1, True, 0], {i, n}, {j, n}];
kms[n_, r_] := Table[r^Abs[i - j], {i, n}, {j, n}];
bhank[n_] := Table[BernoulliB[i + j - 2], {i, n}, {j, n}];
fhank[n_] := Table[(i + j - 2)!, {i, n}, {j, n}];
lauchliA[n_, eps_] := Prepend[eps IdentityMatrix[n], Table[1, {n}]];
lauchliB[n_] := Prepend[Range[n], 1];
(* Dorr: h = 1/(n+1), m = Floor[(n+1)/2], t_i = 1/2 - i h; for i <= m the
   sub-diagonal c_i = -theta/h^2 and super-diagonal e_i = c_i - t_i/h, for
   i > m e_i = -theta/h^2 and c_i = e_i + t_i/h; diagonal -(c_i + e_i). *)
dorr[n_, th_] := Module[{h = 1/(n + 1), m = Floor[(n + 1)/2], cc, ee},
  cc = Table[If[i <= m, -th/h^2, -th/h^2 + (1/2 - i h)/h], {i, n}];
  ee = Table[If[i <= m, -th/h^2 - (1/2 - i h)/h, -th/h^2], {i, n}];
  Table[Which[i == j, -(cc[[i]] + ee[[i]]), j == i - 1, cc[[i]], j == i + 1, ee[[i]], True, 0],
        {i, n}, {j, n}]];

(* Numeric matrices are built OUTSIDE the timed expression: the rows time the
   factorisation, not N of an exact matrix (which a system may cache). *)

(* ---- c1.1 Vandermonde --------------------------------------------------- *)
run["c1.1 Det exact n=30", Det[vand[30]], exactChk];
m = N[vand[8]];
run["c1.1 SingularValueList machine n=8", kappa[m], sig[#, 5] &, 3];
m = N[vand[60], 300];
run["c1.1 SingularValueList wp300 n=60", kappa[m], sig[#, 20] &];

(* ---- c1.2 Gaussian Toeplitz -------------------------------------------- *)
run["c1.2 Det exact n=40", Det[gtoep[40, 9/10]], ratChk];
m = N[gtoep[30, 9/10]];
run["c1.2 Eigenvalues machine n=30", lmin[m], sig[#, 5] &, 3];
m = N[gtoep[120, 19/20], 60];
run["c1.2 Eigenvalues wp60 n=120", lmin[m], sig[#, 20] &];

(* ---- c1.3 Moler --------------------------------------------------------- *)
run["c1.3 Inverse exact n=30", Inverse[moler[30]][[1, 1]], exactChk];
m = N[moler[12]];
run["c1.3 SingularValueList machine n=12", kappa[m], sig[#, 5] &, 3];
m = N[moler[150], 250];
run["c1.3 SingularValueList wp250 n=150", kappa[m], sig[#, 20] &];

(* ---- c1.4 grounded weighted path Laplacian ----------------------------- *)
run["c1.4 LinearSolve exact n=60",
    (Function[x, If[VectorQ[x, IntegerQ], {x[[60]], x[[30]]}, "unevaluated"]])[
      LinearSolve[lap[60], Table[If[i == 60, 1, 0], {i, 60}]]], # &];
m = N[lap[20]];
run["c1.4 Eigenvalues machine n=20", symKappa[m], sig[#, 6] &, 3];
m = N[lap[60], 60];
run["c1.4 Eigenvalues wp60 n=60", symKappa[m], sig[#, 20] &];

(* ---- c1.5 near-singular circulant -------------------------------------- *)
m = N[circ[32, 10^-6]];
run["c1.5 SingularValueList machine n=32", kappa[m], sig[#, 6] &, 3];
m = N[circ[48, 10^-30], 60];
run["c1.5 SingularValueList wp60 n=48", kappa[m], sig[#, 20] &];

(* ---- c1.6 Kac-Murdock-Szego -------------------------------------------- *)
run["c1.6 Det exact n=50", Det[kms[50, 999/1000]], ratChk];
m = N[kms[50, 999/1000]];
run["c1.6 Eigenvalues machine n=50", lmin[m], sig[#, 7] &, 3];
m = N[kms[50, 1 - 10^-12], 60];
run["c1.6 Eigenvalues wp60 n=50", lmin[m], sig[#, 20] &];

(* ---- c1.7 Bernoulli Hankel --------------------------------------------- *)
run["c1.7 Det exact n=14", Det[bhank[14]], ratChk];
m = N[bhank[8]];
run["c1.7 SingularValueList machine n=8", kappa[m], sig[#, 8] &, 3];
m = N[bhank[20], 60];
run["c1.7 SingularValueList wp60 n=20", kappa[m], sig[#, 20] &];

(* ---- c1.8 Lauchli least squares ---------------------------------------- *)
run["c1.8 LeastSquares exact n=20", LeastSquares[lauchliA[20, 10^-10], lauchliB[20]][[20]], ratChk];
{m, rhs} = {N[lauchliA[20, 10^-6]], N[lauchliB[20]]};
run["c1.8 LeastSquares machine n=20", LeastSquares[m, rhs][[20]], sig[#, 6] &, 3];
{m, rhs} = {N[lauchliA[20, 10^-30], 80], N[lauchliB[20], 80]};
run["c1.8 LeastSquares wp80 n=20", LeastSquares[m, rhs][[20]], sig[#, 20] &];

(* ---- c1.9 factorial Hankel --------------------------------------------- *)
run["c1.9 Det exact n=30", Det[fhank[30]], exactChk];
m = N[fhank[6]];
run["c1.9 SingularValueList machine n=6", kappa[m], sig[#, 5] &, 3];
m = N[fhank[60], 400];
run["c1.9 SingularValueList wp400 n=60", kappa[m], sig[#, 20] &];

(* ---- c1.10 Dorr --------------------------------------------------------- *)
run["c1.10 LinearSolve exact n=50", LinearSolve[dorr[50, 1/100], Table[1, {50}]][[25]], sig[#, 20] &];
m = N[dorr[60, 1/200]];
run["c1.10 LinearSolve machine n=60", LinearSolve[m, Table[1., {60}]][[30]], sig[#, 5] &, 3];
m = N[dorr[100, 1/1000], 60];
run["c1.10 LinearSolve wp60 n=100", LinearSolve[m, Table[N[1, 60], {100}]][[50]], sig[#, 20] &];
