(* bench-timeout: 2400 *)
(* ==========================================================================
   Experiment 164 -- ISSAC '97 holdout C (SEALED), family c9:
   eigenvalues of operators and of their discretisations (counterpart of
   Problem 9)
   ==========================================================================

       cd benchmarks/164-issac97c-operator-eigenvalues
       ../../Mathilda -file operator_eigenvalues.m

   Ten eigenproblems: three NON-NORMAL matrices (a tridiagonal Toeplitz
   matrix similar to a symmetric one only through a diagonal scaling of
   condition 1e236, a banded Toeplitz matrix with symbol z + z^-2/4, and a
   Jordan block closed by eps = 3e-60), two Sturm-Liouville problems with a
   SINGULAR endpoint (a Bessel operator of order 1/3 at x = 0, and the
   linear potential |x| on the whole line), a regular radial problem on an
   annulus, two Fredholm operators (Love's kernel 1/(Pi (1 + (x-y)^2)) on
   [-1, 1] by Gauss-Legendre Nystrom with Golub-Welsch nodes, and e^(xy) on
   [0, 1] by its exact degenerate-kernel expansion), the PT-symmetric
   non-Hermitian oscillator p^2 + i x^3 in a harmonic-oscillator basis
   (complex symmetric), and the algebraic connectivity of a lollipop graph.

   Machine rows: the discretisation (or shooting) in doubles; wp rows: the
   same at 40-300 digits, larger where the discretisation must converge
   further.  Shooting rows (c9.2, c9.6, c9.10) are NDSolve inside FindRoot,
   checked to 6 digits (machine) / 15 digits (WorkingPrecision -> 40).

   TRUE ANSWERS (A and B independent, agreeing to >= 25 digits unless noted):
   c9.1  T = tridiag(eps, 1, 1), n x n: lambda_max = 1 + 2 Sqrt[eps] Cos[Pi/(n+1)]
         (closed form; diagonal similarity to a symmetric matrix).  machine
         n=6, eps=1/100: 1.18019377358048382524722046390; wp300 n=60,
         eps=1e-8: 1.00019973481797696610152115435.
   c9.2  -(x u')' + (1/9) u/x = lambda x u on (0, 1), u(1) = 0, u regular at 0
         (Frobenius start u = x^(1/3) (1 - 3 lambda x^2/16)): lambda =
         BesselJZero[1/3, 1]^2 = 8.42500692949919857451071877294 (A mpmath
         besseljzero; B PARI solve of besselj).
   c9.3  Love: lambda_max of (1/Pi) int_-1^1 u(y)/(1 + (x-y)^2) dy
         = 0.453572864122606855682199755895 (A mpmath Nystrom, Golub-Welsch
         Gauss-Legendre at n = 40/50/60 agree to 1e-34; machine row n = 20 is
         converged to 5e-18).
   c9.4  H = p^2 + i x^3, ground state, oscillator basis (frequency 3, the x^3
         matrix formed in a basis 4 larger and truncated): E0 =
         1.15626707198811325616 (A mpmath eig at N = 80/100/120; B the
         published 1.1562670719881132, Bender & Boettcher 1998).  machine
         N = 40 (error 2.7e-11), wp60 N = 100.
   c9.5  Lollipop: K_10 with a 10-vertex path attached (20 vertices), the
         second-smallest Laplacian eigenvalue = 0.0379599847847906607432036502309
         (A mpmath eigsy; B PARI polroots of the exact characteristic polynomial).
   c9.6  -u'' + |x| u = E u on the line, ground state: E0 = -AiryAiPrime zero
         = 1.01879297164747108901732478340 (A mpmath; B PARI).  Shooting on
         [0, L] with u'(0) = 0 (L = 12 machine, 16 wp: the truncation error is
         e^(-(4/3)(L-E)^(3/2)) < 1e-21).
   c9.7  e^(xy) on [0, 1]: e^(xy) = sum x^j y^j/j!, so the operator acts on
         span{x^j} by M_jk = 1/(j! (j+k+1)), j, k = 0..N-1.  lambda_max
         = 1.35303016474573485330920174962 (A mpmath eig, N = 40; B Nystrom
         Gauss-Legendre n = 30, same to 34 digits).  machine N = 20, wp60 N = 40.
   c9.8  T_(i,i-1) = 1, T_(i,i+2) = 1/4 (non-normal banded Toeplitz):
         spectral radius.  machine n = 12: 1.12491323229113313151398617253;
         wp200 n = 40: 1.18367763240231769750996446175 (A mpmath eig at 60/300
         digits; B PARI polroots of the exact characteristic polynomial).
   c9.9  Jordan block J_n(0) + eps e_n e_1^T: the characteristic polynomial is
         lambda^n - eps, so the spectral radius is eps^(1/n) (closed form).
         machine n = 8, eps = 2e-8: 0.109050773266525765920701065576;
         wp200 n = 30, eps = 3e-60: 0.0103729919730295071056672534457.
   c9.10 -(r u')' = lambda r u on (1, 2), u(1) = u(2) = 0: lambda = k^2 with
         J0(k) Y0(2k) = J0(2k) Y0(k): 9.75332212475071491068952390873 (A mpmath
         findroot; B PARI solve with besselj/besseln).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Eigenvalues", "Eigensystem", "NDSolve", "FindRoot"}];

lamChk[d_][r_] := Module[{v = ruleVal[r, lam]}, If[NumericQ[v], sig[Re[v], d], "failed"]];

(* ---- c9.1 non-normal tridiagonal Toeplitz ------------------------------ *)
mat1[n_, e_] := Table[Which[i == j, 1, j == i + 1, 1, i == j + 1, e, True, 0], {i, n}, {j, n}];
m = N[mat1[6, 1/100]];
run["c9.1 Eigenvalues machine n=6", Max[Re[Eigenvalues[m]]], sig[#, 8] &, 3];
m = N[mat1[60, 10^-8], 300];
run["c9.1 Eigenvalues wp300 n=60", Max[Re[Eigenvalues[m]]], sig[#, 20] &];

(* ---- c9.2 Bessel operator of order 1/3, shooting from the singular end -- *)
u2[l_?NumericQ, x0_, wp_] := Module[{s},
  s = NDSolve[{u''[x] == -u'[x]/x + u[x]/(9 x^2) - l u[x],
               u[x0] == x0^(1/3) (1 - 3 l x0^2/16), u'[x0] == x0^(-2/3)/3 - 7 l x0^(4/3)/16},
              u, {x, x0, 1}, WorkingPrecision -> wp, MaxSteps -> 100000];
  u[1] /. First[s]];
run["c9.2 shooting machine", FindRoot[u2[lam, 10^-4, MachinePrecision], {lam, 8, 9}], lamChk[6]];
run["c9.2 shooting wp40", FindRoot[u2[lam, 10^-6, 40], {lam, 8, 9}, WorkingPrecision -> 40], lamChk[15]];

(* ---- c9.3 Love's kernel by Gauss-Legendre Nystrom ----------------------- *)
gl[n_, prec_] := Module[{es = Eigensystem[N[Table[Which[j == i + 1, i/Sqrt[4 i^2 - 1],
                                                      i == j + 1, j/Sqrt[4 j^2 - 1], True, 0],
                                                {i, n}, {j, n}], prec]]},
  {es[[1]], Map[2 #[[1]]^2/(# . #) &, es[[2]]]}];
love[n_, prec_] := Module[{xw = gl[n, prec], xs, ws},
  xs = xw[[1]]; ws = xw[[2]];
  Max[Eigenvalues[Table[Sqrt[ws[[i]] ws[[j]]]/(N[Pi, prec] (1 + (xs[[i]] - xs[[j]])^2)), {i, n}, {j, n}]]]];
run["c9.3 Nystrom machine n=20", love[20, MachinePrecision], sig[#, 8] &, 3];
run["c9.3 Nystrom wp50 n=50", love[50, 50], sig[#, 20] &];

(* ---- c9.4 PT-symmetric p^2 + i x^3 in the oscillator basis -------------- *)
ptE[nb_, w_, prec_] := Module[{X, P2, H},
  X = N[Table[Which[j == i + 1, Sqrt[i/(2 w)], i == j + 1, Sqrt[j/(2 w)], True, 0], {i, nb + 4}, {j, nb + 4}], prec];
  P2 = N[Table[Which[i == j, w (2 i - 1)/2, j == i + 2, -w/2 Sqrt[i (i + 1)],
                     i == j + 2, -w/2 Sqrt[j (j + 1)], True, 0], {i, nb + 4}, {j, nb + 4}], prec];
  H = P2 + I X . X . X;
  Min[Re[Eigenvalues[Table[H[[i, j]], {i, nb}, {j, nb}]]]]];
run["c9.4 oscillator basis machine N=40", ptE[40, 3, MachinePrecision], sig[#, 8] &, 3];
run["c9.4 oscillator basis wp60 N=100", ptE[100, 3, 60], sig[#, 15] &];

(* ---- c9.5 lollipop graph Laplacian -------------------------------------- *)
lolEdges = Join[Flatten[Table[{i, j}, {i, 10}, {j, i + 1, 10}], 1], Table[{i, i + 1}, {i, 10, 19}]];
lolA = Table[If[MemberQ[lolEdges, {Min[i, j], Max[i, j]}], 1, 0], {i, 20}, {j, 20}];
lolL = DiagonalMatrix[Total /@ lolA] - lolA;
m = N[lolL];
run["c9.5 Eigenvalues machine", Sort[Eigenvalues[m]][[2]], sig[#, 8] &, 3];
m = N[lolL, 40];
run["c9.5 Eigenvalues wp40", Sort[Eigenvalues[m]][[2]], sig[#, 20] &];

(* ---- c9.6 linear potential, shooting on [0, L] -------------------------- *)
u6[e_?NumericQ, L_, wp_] := Module[{s},
  s = NDSolve[{u''[x] == (x - e) u[x], u[0] == 1, u'[0] == 0}, u, {x, 0, L},
              WorkingPrecision -> wp, MaxSteps -> 100000];
  u[L] /. First[s]];
run["c9.6 shooting machine", FindRoot[u6[lam, 12, MachinePrecision], {lam, 1, 105/100}], lamChk[6]];
run["c9.6 shooting wp40", FindRoot[u6[lam, 16, 40], {lam, 1, 105/100}, WorkingPrecision -> 40], lamChk[15]];

(* ---- c9.7 e^(xy): degenerate-kernel matrix ------------------------------ *)
mat7[n_] := Table[1/((j - 1)! (j + k - 1)), {j, n}, {k, n}];
m = N[mat7[20]];
run["c9.7 Eigenvalues machine N=20", Max[Re[Eigenvalues[m]]], sig[#, 8] &, 3];
m = N[mat7[40], 60];
run["c9.7 Eigenvalues wp60 N=40", Max[Re[Eigenvalues[m]]], sig[#, 20] &];

(* ---- c9.8 non-normal banded Toeplitz ------------------------------------ *)
mat8[n_] := Table[Which[i - j == 1, 1, j - i == 2, 1/4, True, 0], {i, n}, {j, n}];
m = N[mat8[12]];
run["c9.8 Eigenvalues machine n=12", Max[Abs[Eigenvalues[m]]], sig[#, 8] &, 3];
m = N[mat8[40], 200];
run["c9.8 Eigenvalues wp200 n=40", Max[Abs[Eigenvalues[m]]], sig[#, 20] &];

(* ---- c9.9 Jordan block closed by eps ------------------------------------ *)
mat9[n_, e_] := Table[Which[j == i + 1, 1, i == n && j == 1, e, True, 0], {i, n}, {j, n}];
m = N[mat9[8, 2 10^-8]];
run["c9.9 Eigenvalues machine n=8", Max[Abs[Eigenvalues[m]]], sig[#, 8] &, 3];
m = N[mat9[30, 3 10^-60], 200];
run["c9.9 Eigenvalues wp200 n=30", Max[Abs[Eigenvalues[m]]], sig[#, 20] &];

(* ---- c9.10 annulus, shooting -------------------------------------------- *)
u10[l_?NumericQ, wp_] := Module[{s},
  s = NDSolve[{u''[x] == -u'[x]/x - l u[x], u[1] == 0, u'[1] == 1}, u, {x, 1, 2},
              WorkingPrecision -> wp, MaxSteps -> 100000];
  u[2] /. First[s]];
run["c9.10 shooting machine", FindRoot[u10[lam, MachinePrecision], {lam, 9, 10}], lamChk[6]];
run["c9.10 shooting wp40", FindRoot[u10[lam, 40], {lam, 9, 10}, WorkingPrecision -> 40], lamChk[15]];
