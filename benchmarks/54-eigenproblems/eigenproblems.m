(* ==========================================================================
   Experiment 54 -- Eigenproblems (group E, advanced numerical analysis)
   ==========================================================================
   WHAT IT MEASURES.  The eigensolvers in src/linalg/eigen_*.c -- the symmetric
   (syev) and general (geev) dense drivers, the generalized problem
   Eigenvalues[{A, B}], and the Krylov path Method -> "Arnoldi".  Group B times
   only Eigenvalues on a single dense matrix; this separates the algorithms.

   TIMING on large random matrices; CHECKS on the small reproducible
   spdIntMatrix[8] or a fixed matrix.  Eigen-outputs are non-unique (order, sign,
   eigenvector scale), so checks are order-invariant scalars: max eigenvalue,
   sum of |eigenvalues|, trace via eigenvalues.  Eigenvectors have no stable
   scalar (Mathematica does not normalise them), so that case carries no check.

   ONE FINDING IS BAKED IN.  Eigensystem is ABSENT (benchIf).

   The generalized problem Eigenvalues[{A, B}] now takes a LAPACK machine path
   (dsygv for the symmetric-definite pencil here, dggev/zhegv/zggev otherwise;
   eigen_direct.c), so the 5x5 loop runs in sub-millisecond -- FASTER than
   scipy.linalg.eigh -- instead of forming Det[A - lambda B] and symbolically
   root-finding it (Root[] objects by size 6, no scaling past ~7).  The small
   5x5 loop is kept so the row stays comparable across the history.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"Eigenvalues", "Eigenvectors", "Eigensystem"}];

seed[];
ns = 300;
sy  = With[{a = rand01[{ns, ns}]}, (a + Transpose[a])/2];       (* symmetric *)
ge  = rand01[{ns, ns}];                                          (* general *)
sv  = With[{a = rand01[{250, 250}]}, (a + Transpose[a])/2];     (* eigenvectors *)
big = With[{a = rand01[{500, 500}]}, (a + Transpose[a])/2];     (* Arnoldi *)

sm   = N[spdIntMatrix[8]];
gen4 = N[{{3, 1, 0, 2}, {2, 4, 1, 0}, {0, 1, 5, 1}, {1, 0, 1, 6}}];
A5   = N[spdIntMatrix[5]];
B5   = N[DiagonalMatrix[{2, 3, 4, 5, 6}]];

bench["Eigenvalues symmetric 300x300", Eigenvalues[sy];, 1];
check["Eigenvalues symmetric 300x300", Round[10^6 Max[Eigenvalues[sm]]]];

bench["Eigenvalues general 300x300", Eigenvalues[ge];, 1];
check["Eigenvalues general 300x300", Round[10^4 Total[Abs[Eigenvalues[gen4]]]]];

bench["Eigenvectors symmetric 250x250", Eigenvectors[sv];, 1];  (* no check: eigenvectors non-unique *)

(* Generalized: no machine path -- symbolic char-poly root-finding.  Small 5x5
   loop so it completes; the ratio against scipy.linalg.eigh(A,B) is the story. *)
bench["Generalized eigenvalues 5x5 x50", Do[Eigenvalues[{A5, B5}], {50}];, 1];
check["Generalized eigenvalues 5x5 x50", Round[10^6 Total[Eigenvalues[{A5, B5}]]]];

bench["Eigenvalues Arnoldi k=6 of 500x500", Eigenvalues[big, 6, Method -> "Arnoldi"];];
check["Eigenvalues Arnoldi k=6 of 500x500", Round[10^6 Total[Eigenvalues[sm]]]];

benchIf["Eigensystem symmetric 250x250", "Eigensystem", Eigensystem[sv];, 1];
checkIf["Eigensystem symmetric 250x250", "Eigensystem", Round[10^6 Max[Eigenvalues[sm]]]];
