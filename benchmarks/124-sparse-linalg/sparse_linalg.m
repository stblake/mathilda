(* ==========================================================================
   Experiment 124 -- sparse linear algebra at scale
   ==========================================================================
   WHAT IT MEASURES.  The native sparse solvers of src/sparse/ (v0.6180-
   v0.6188) on the problems they exist for, one size class up from
   120-sparse-arrays:

       poisson2d   5-point Laplacian, 500^2 = 2.5 x 10^5 unknowns, direct
       poisson3d   7-point Laplacian, 30^3 = 2.7 x 10^4 direct (3-D fill),
                   100^3 = 10^6 unknowns by unpreconditioned CG
       convdiff    nonsymmetric convection-diffusion on 300^2: UMFPACK,
                   GMRES + ILUT, BiCGSTAB + ILU0
       helmholtz   complex-shifted indefinite Laplacian on 300^2, direct
       lsf         one factorisation, 100 right-hand sides (LinearSolveFunction)
       heat        one implicit-free heat step MatrixExp[-L, u0] on 500^2
       eigen       6 largest / 6 smallest eigenvalues and Eigensystem of a
                   200^2 Laplacian plus a separating potential (ARPACK) --
                   last in the file, see below
       exact       small exact integer matrices (10 x 10 .. 40 x 40, a third
                   of the entries nonzero): Det, Inverse, LinearSolve,
                   MatrixRank, NullSpace, CharacteristicPolynomial

   Matrices are built from Band and KroneckerProduct, so every system builds
   the same operator; right-hand sides are (k mod 7) - 3.  The exact matrices
   come from a Lehmer step s = (48271 (n (i - 1) + j) + seed) mod (2^31 - 1):
   entry (s mod 19) - 9 where s mod 3 = 0, else 0, plus 10 on the diagonal.
   A rank-r matrix is exM[n, 7] . G[n x r] . G[r x n] of two more of them.

   CHECKS.  Numeric solves check Round[10 Norm[x]] or Round[Norm[x]], the
   scale picked per case so the value sits well clear of a rounding boundary
   (a 10^-8 Krylov residual moves Total[x], which cancels, but not the norm);
   the heat step checks Round[1000 Norm]; eigenvalues Round[10^6 Total].  The exact cases check the
   exact answer: the determinant, the sum of the inverse's entries and of the
   solution (a Rational), the rank, the nullity, and the characteristic
   polynomial at x = 3 (n is even, so Det[m - x I] and det(x I - m) agree).

   Sizes are held where the scipy column still fits the file budget: its
   direct solver is SuperLU, which costs ~8x Mathilda's CHOLMOD on the 500^2
   grid and ~100x on 3-D fill (124 s at 40^3, against 1.5 s), so the 3-D
   direct row is 30^3 and the 10^6-unknown row is the 3-D CG.
   ========================================================================== *)

Get["../harness.m"];

require[{"SparseArray", "LinearSolve", "LinearSolveFunction", "MatrixExp",
         "Eigenvalues", "Eigensystem", "Det", "Inverse", "MatrixRank",
         "NullSpace", "CharacteristicPolynomial"}];

tri[k_] := SparseArray[{Band[{1, 1}] -> 2., Band[{1, 2}] -> -1., Band[{2, 1}] -> -1.}, {k, k}];
eye[k_] := IdentityMatrix[k, SparseArray];
lap2[k_] := KroneckerProduct[tri[k], eye[k]] + KroneckerProduct[eye[k], tri[k]];
lap3[k_] := KroneckerProduct[tri[k], eye[k^2]] + KroneckerProduct[eye[k], tri[k], eye[k]] +
            KroneckerProduct[eye[k^2], tri[k]];
rhs[n_] := N[Mod[Range[n], 7]] - 3.;
nrm10[x_] := Round[10 Norm[x]];
krylov[a_, b_, meth_, pre_, tol_] := LinearSolve[a, b,
  Method -> {"Krylov", Method -> meth, "Preconditioner" -> pre, Tolerance -> tol}];

(* ---- Poisson ------------------------------------------------------------- *)

L500 = lap2[500]; b500 = rhs[500^2];
bench["poisson2d 500^2: LinearSolve direct", x = LinearSolve[L500, b500], 1];
check["poisson2d 500^2: LinearSolve direct", nrm10[x]];

L30 = lap3[30]; b30 = rhs[30^3];
benchOnce["poisson3d 30^3: LinearSolve direct", x = LinearSolve[L30, b30]];
check["poisson3d 30^3: LinearSolve direct", nrm10[x]];

L100 = lap3[100]; b100 = rhs[100^3];
benchOnce["poisson3d 100^3 (10^6): CG 1e-6",
          x = krylov[L100, b100, "ConjugateGradient", None, 10^-6]];
check["poisson3d 100^3 (10^6): CG 1e-6", Round[Norm[x]]];

(* ---- convection-diffusion (nonsymmetric) ---------------------------------- *)

skew = SparseArray[{Band[{1, 2}] -> 1., Band[{2, 1}] -> -1.}, {300, 300}];
C300 = lap2[300] + 0.4 KroneckerProduct[eye[300], skew] + 0.2 KroneckerProduct[skew, eye[300]];
b300 = rhs[300^2];
bench["convdiff 300^2: LinearSolve Multifrontal", x = LinearSolve[C300, b300, Method -> "Multifrontal"], 1];
check["convdiff 300^2: LinearSolve Multifrontal", Round[Norm[x]]];
bench["convdiff 300^2: GMRES + ILUT 1e-8", x = krylov[C300, b300, "GMRES", "ILUT", 10^-8], 1];
check["convdiff 300^2: GMRES + ILUT 1e-8", Round[Norm[x]]];
bench["convdiff 300^2: BiCGSTAB + ILU0 1e-8", x = krylov[C300, b300, "BiCGSTAB", "ILU0", 10^-8], 1];
check["convdiff 300^2: BiCGSTAB + ILU0 1e-8", Round[Norm[x]]];

(* ---- complex Helmholtz-like: shifted into the spectrum, lightly damped ---- *)

H300 = lap2[300] + (-0.3 + 0.05 I) eye[300^2];
bench["helmholtz 300^2 complex: LinearSolve direct", x = LinearSolve[H300, b300], 1];
check["helmholtz 300^2 complex: LinearSolve direct", nrm10[x]];

(* ---- many right-hand sides: factor once (inside the timing), 100 solves --- *)

L300 = lap2[300];
bench["lsf 300^2: LinearSolveFunction, factor + 100 solves",
      (f = LinearSolve[L300]; Do[x = f[b300 + k], {k, 100}]), 1];
check["lsf 300^2: LinearSolveFunction, factor + 100 solves", Round[Norm[f[b300 + 100.]]]];

(* ---- heat step: the action of exp(-L) on u0, no dense exponential -------- *)

u500 = Sin[N[Range[500^2]]/100.];
bench["heat 500^2: MatrixExp[-L, u0] (action)", y = MatrixExp[-L500, u500]];
check["heat 500^2: MatrixExp[-L, u0] (action)", Round[1000 Norm[y]]];

(* ---- exact small matrices ------------------------------------------------- *)

exE[n_, sd_, i_, j_] := With[{s = Mod[48271 (n (i - 1) + j) + sd, 2147483647]},
  If[Mod[s, 3] == 0, Mod[s, 19] - 9, 0] + If[i == j, 10, 0]];
exM[n_, sd_] := Table[exE[n, sd, i, j], {i, n}, {j, n}];
exR[n_, r_] := exM[n, 7] . Table[exE[n, 11, i, j], {i, n}, {j, r}] .
               Table[exE[n, 13, i, j], {i, r}, {j, n}];

m40 = exM[40, 1]; m30 = exM[30, 2]; m20 = exM[20, 3]; m10 = exM[10, 4];
r40 = exR[40, 30]; r30 = exR[30, 24];

bench["exact 40x40: Det", d = Det[m40]];
check["exact 40x40: Det", d];
bench["exact 30x30: Inverse", iv = Inverse[m30]];
check["exact 30x30: Inverse", Total[iv, 2]];
bench["exact 40x40: LinearSolve", x = LinearSolve[m40, Range[40]]];
check["exact 40x40: LinearSolve", Total[x]];
bench["exact 40x40: MatrixRank (rank 30)", rk = MatrixRank[r40]];
check["exact 40x40: MatrixRank (rank 30)", rk];
bench["exact 30x30: NullSpace (nullity 6)", ns = NullSpace[r30]];
check["exact 30x30: NullSpace (nullity 6)", Length[ns]];
bench["exact 20x20: CharacteristicPolynomial", cp = CharacteristicPolynomial[m20, xx]];
check["exact 20x20: CharacteristicPolynomial", cp /. xx -> 3];
bench["exact 10x10: Inverse x100", Do[iv = Inverse[m10], {100}]];
check["exact 10x10: Inverse x100", Total[iv, 2]];

(* ---- eigenvalues with a count (ARPACK) -----------------------------------
   Last on purpose: the sparse eigensolvers are still being built, and a
   timeout here loses only these rows.  H = L + V with the potential
   V(i, j) = 0.001 (i^2 + 2 j^2): the bare Laplacian's extreme eigenvalues sit
   ~3e-4 apart (relative) in pairs, and at that gap the reference's default
   Arnoldi tolerance returns visibly unconverged values (see 120); the
   potential separates both ends of the spectrum so every system converges to
   the same six. *)

H200 = lap2[200] + SparseArray[Band[{1, 1}] -> 0.001 Flatten[Table[i^2 + 2 j^2, {i, 200}, {j, 200}]],
                               {200^2, 200^2}];
bench["eigen 200^2 L + V: Eigenvalues[H, 6] (largest)", ev = Eigenvalues[H200, 6], 1];
check["eigen 200^2 L + V: Eigenvalues[H, 6] (largest)", Round[10^6 Total[ev]]];
bench["eigen 200^2 L + V: Eigenvalues[H, -6] (smallest)", ev = Eigenvalues[H200, -6], 1];
check["eigen 200^2 L + V: Eigenvalues[H, -6] (smallest)", Round[10^6 Total[ev]]];
bench["eigen 200^2 L + V: Eigensystem[H, -6]", es = Eigensystem[H200, -6], 1];
check["eigen 200^2 L + V: Eigensystem[H, -6]",
      {Round[10^6 Total[First[es]]], Round[1000 Total[Flatten[Abs[Last[es]]^2]]]}];
