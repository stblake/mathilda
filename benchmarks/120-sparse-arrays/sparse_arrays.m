(* ==========================================================================
   Experiment 120 -- SparseArray construction, arithmetic and element-wise
   functions against scipy.sparse
   ==========================================================================
   WHAT IT MEASURES.  The sparse-array substrate a sparse linear-algebra or
   graph workload sits on, on a 10^5 x 10^5 matrix with 10^6 explicit
   entries (density 10^-4):

       build      coordinate list -> compressed rows   (SparseArray[pos -> v])
       arith      A + B, A B (element-wise), 2.5 A, A^2, A + A
       functions  Sin[A], Sqrt[A]                      (sparsity-preserving)
       access     one row, one element, the explicit entries
       dense      a 2000 x 2000 dense <-> sparse round trip

   Positions are a Lehmer LCG, i_k = (48271 k mod (2^31 - 1)) mod n + 1 and
   j_k = (69621 k mod (2^31 - 1)) mod n + 1 (B shifts j by 7 before the
   reduction for odd k only, so A and B share half their positions), so
   every system builds the same matrices; no position repeats within a
   matrix, so first-rule-wins and sum-duplicates agree.  Values are
   (k mod 97)/10 + 0.1 > 0.

   CHECK.  Round[1000 x (sum of the explicit values)] of each result -- an
   Integer, so the systems compare it exactly.
   ========================================================================== *)

Get["../harness.m"];

require[{"SparseArray", "ArrayRules", "Normal"}];

$n = 100000; $nnz = 1000000;
kk = Range[$nnz];
ii = Mod[Mod[48271 kk, 2147483647], $n] + 1;
jj = Mod[Mod[69621 kk, 2147483647], $n] + 1;
jj2 = Mod[Mod[69621 kk + 7 Mod[kk, 2], 2147483647], $n] + 1;
vv = N[Mod[kk, 97]]/10. + 0.1;
posA = Transpose[{ii, jj}];
posB = Transpose[{ii, jj2}];

digest[s_] := Round[1000 Total[s["NonzeroValues"]]];

(* ---- build -------------------------------------------------------------- *)

bench["build: 10^6 entries, 10^5 x 10^5 (coordinates -> rows)",
      A = SparseArray[posA -> vv, {$n, $n}]];
check["build: 10^6 entries, 10^5 x 10^5 (coordinates -> rows)", digest[A]];
B = SparseArray[posB -> vv, {$n, $n}];

(* ---- arithmetic --------------------------------------------------------- *)

bench["arith: A + B", r = A + B];
check["arith: A + B", digest[r]];
bench["arith: A B (element-wise)", r = A B];
check["arith: A B (element-wise)", digest[r]];
bench["arith: 2.5 A", r = 2.5 A];
check["arith: 2.5 A", digest[r]];
bench["arith: A^2 (element-wise)", r = A^2];
check["arith: A^2 (element-wise)", digest[r]];
bench["arith: A + A", r = A + A];
check["arith: A + A", digest[r]];

(* ---- element-wise functions --------------------------------------------- *)

bench["fn: Sin[A]", r = Sin[A]];
check["fn: Sin[A]", digest[r]];
bench["fn: Sqrt[A]", r = Sqrt[A]];
check["fn: Sqrt[A]", digest[r]];

(* ---- access ------------------------------------------------------------- *)

bench["access: row 500 x100", Do[r = A[[500]], {100}]];
check["access: row 500 x100", Round[1000 Total[Normal[A[[500]]]]]];
bench["access: element (500, 600) x1000", Do[r = A[[500, 600]], {1000}]];
bench["access: explicit positions and values", (p = A["NonzeroPositions"]; v = A["NonzeroValues"])];
check["access: explicit positions and values", Total[p, 2]];

(* ---- dense <-> sparse --------------------------------------------------- *)

dense = Table[If[Mod[i + 3 j, 97] == 0, 1.5 + Mod[i j, 7], 0.], {i, 2000}, {j, 2000}];
bench["dense: 2000 x 2000 dense -> sparse", Dm = SparseArray[dense]];
check["dense: 2000 x 2000 dense -> sparse", digest[Dm]];
bench["dense: 2000 x 2000 sparse -> dense", r = Normal[Dm]];
check["dense: 2000 x 2000 sparse -> dense", Round[1000 Total[r, 2]]];

(* ---- products ------------------------------------------------------------ *)

x = N[Mod[Range[$n], 13]]/7.;
bench["dot: A . x (SpMV) x100", Do[r = A . x, {100}]];
check["dot: A . x (SpMV) x100", Round[1000 Total[A . x]]];
bench["dot: x . A x10", Do[r = x . A, {10}]];
check["dot: x . A x10", Round[1000 Total[x . A]]];
bench["dot: A . A (SpGEMM)", r = A . A];
check["dot: A . A (SpGEMM)", digest[r]];
bench["struct: Transpose[A]", At = Transpose[A]];
check["struct: Transpose[A]", digest[At]];
bench["dot: A . Transpose[A] (SpGEMM)", r = A . At];
check["dot: A . Transpose[A] (SpGEMM)", digest[r]];

(* ---- reductions ----------------------------------------------------------- *)

bench["reduce: Total[A] (column sums)", r = Total[A]];
check["reduce: Total[A] (column sums)", Round[1000 Total[r]]];
bench["reduce: Total[A, {2}] (row sums)", r = Total[A, {2}]];
check["reduce: Total[A, {2}] (row sums)", Round[1000 Total[r]]];
bench["reduce: Total[A, 2]", r = Total[A, 2]];
check["reduce: Total[A, 2]", Round[1000 r]];
bench["reduce: Max[A]", r = Max[A]];
check["reduce: Max[A]", Round[1000 r]];
bench["reduce: Norm[A, Frobenius]", r = Norm[A, "Frobenius"]];
check["reduce: Norm[A, Frobenius]", Round[1000 r]];

(* ---- assembly ------------------------------------------------------------- *)

bench["struct: Join[A, B] (vstack)", r = Join[A, B]];
check["struct: Join[A, B] (vstack)", digest[r]];
bench["struct: ArrayFlatten 2x2 blocks", r = ArrayFlatten[{{A, B}, {B, A}}]];
check["struct: ArrayFlatten 2x2 blocks", digest[r]];
bench["struct: IdentityMatrix[10^6, SparseArray]", r = IdentityMatrix[10^6, SparseArray]];
check["struct: IdentityMatrix[10^6, SparseArray]", Round[Total[r["NonzeroValues"]]]];

(* ---- 2-D Poisson: 5-point Laplacian on a 300 x 300 grid ------------------- *)

$m = 300;
bench["poisson: build L (Band, KroneckerProduct)",
      (T1 = SparseArray[{Band[{1, 1}] -> 2., Band[{1, 2}] -> -1., Band[{2, 1}] -> -1.}, {$m, $m}];
       I1 = IdentityMatrix[$m, SparseArray];
       L = KroneckerProduct[T1, I1] + KroneckerProduct[I1, T1])];
check["poisson: build L (Band, KroneckerProduct)", digest[L]];
bb = N[Mod[Range[$m^2], 7]] - 3.;
bench["poisson: LinearSolve[L, b] direct", sol = LinearSolve[L, bb]];
check["poisson: LinearSolve[L, b] direct", Round[1000 Total[sol]]];
bench["poisson: LinearSolve[L, b] CG 1e-10",
      solcg = LinearSolve[L, bb, Method -> {"Krylov", Method -> "ConjugateGradient",
                                            "Preconditioner" -> None, Tolerance -> 10^-10}]];
check["poisson: LinearSolve[L, b] CG 1e-10", Round[Total[solcg]]];
(* ---- native sparse solvers (v0.6180-v0.6188) on the same L ---------------

   The direct methods by name, the Krylov methods with each incomplete-LU
   preconditioner, a factor reused over 20 right-hand sides
   (LinearSolveFunction), Det / Inverse of a 900 x 900 Laplacian, the action
   MatrixExp[s, v] (no dense exponential formed), MatrixExp of a sparse
   1-D Laplacian, and a complex-shifted (Helmholtz-like) solve.
   Krylov checks are Round[10 Norm[x]]: Total[x] cancels to -56.6 over
   90000 entries, so a 10^-10 residual still moves it in the first decimal,
   while the norm (330.55) is stable to 10^-4. *)

nrm10[x_] := Round[10 Norm[x]];
krylov[meth_, pre_] := LinearSolve[L, bb,
  Method -> {"Krylov", Method -> meth, "Preconditioner" -> pre, Tolerance -> 10^-10}];

bench["poisson: LinearSolve Multifrontal", sol = LinearSolve[L, bb, Method -> "Multifrontal"]];
check["poisson: LinearSolve Multifrontal", Round[1000 Total[sol]]];
bench["poisson: LinearSolve Cholesky", sol = LinearSolve[L, bb, Method -> "Cholesky"]];
check["poisson: LinearSolve Cholesky", Round[1000 Total[sol]]];

bench["poisson: Krylov CG + ILU0 1e-10", x = krylov["ConjugateGradient", "ILU0"], 1];
check["poisson: Krylov CG + ILU0 1e-10", nrm10[x]];
bench["poisson: Krylov BiCGSTAB + ILU0 1e-10", x = krylov["BiCGSTAB", "ILU0"], 1];
check["poisson: Krylov BiCGSTAB + ILU0 1e-10", nrm10[x]];
bench["poisson: Krylov BiCGSTAB + ILUT 1e-10", x = krylov["BiCGSTAB", "ILUT"], 1];
check["poisson: Krylov BiCGSTAB + ILUT 1e-10", nrm10[x]];
bench["poisson: Krylov GMRES + ILUT 1e-10", x = krylov["GMRES", "ILUT"], 1];
check["poisson: Krylov GMRES + ILUT 1e-10", nrm10[x]];
bench["poisson: Krylov GMRES + ILUTP 1e-10", x = krylov["GMRES", "ILUTP"], 1];
check["poisson: Krylov GMRES + ILUTP 1e-10", nrm10[x]];

(* Factor once, then 20 distinct right-hand sides; the factorisation is
   inside the timing, so a cached factor cannot flatter a rerun. *)
bench["poisson: LinearSolveFunction, factor + 20 solves",
      (lsf = LinearSolve[L]; Do[x = lsf[bb + k], {k, 20}])];
check["poisson: LinearSolveFunction, factor + 20 solves", Round[1000 Total[lsf[bb + 20.]]]];

u0 = Sin[N[Range[$m^2]]/100.];
bench["poisson: MatrixExp[-0.5 L, u] (action)", y = MatrixExp[-0.5 L, u0]];
check["poisson: MatrixExp[-0.5 L, u] (action)", Round[1000 Total[y]]];

Lz = L + 0.5 I IdentityMatrix[$m^2, SparseArray];
bench["poisson: complex shift LinearSolve[L + 0.5 I, b]", z = LinearSolve[Lz, bb]];
check["poisson: complex shift LinearSolve[L + 0.5 I, b]",
      {Round[1000 Re[Total[z]]], Round[1000 Im[Total[z]]]}];

T30 = SparseArray[{Band[{1, 1}] -> 2., Band[{1, 2}] -> -1., Band[{2, 1}] -> -1.}, {30, 30}];
I30 = IdentityMatrix[30, SparseArray];
L30 = KroneckerProduct[T30, I30] + KroneckerProduct[I30, T30];
bench["sparse: Det of 900x900 Laplacian / 4", d = Det[L30/4.]];
check["sparse: Det of 900x900 Laplacian / 4", Round[10^6 Log[d]]];
bench["sparse: Inverse of 900x900 Laplacian", inv = Inverse[L30]];
check["sparse: Inverse of 900x900 Laplacian", Round[1000 Total[inv, 2]]];

T1k = SparseArray[{Band[{1, 1}] -> 2., Band[{1, 2}] -> -1., Band[{2, 1}] -> -1.}, {1000, 1000}];
bench["sparse: MatrixExp[-0.1 T] 1-D Laplacian 1000x1000", e = MatrixExp[-0.1 T1k]];
check["sparse: MatrixExp[-0.1 T] 1-D Laplacian 1000x1000", Round[1000 Total[e, 2]]];

bench["poisson: 6 largest eigenvalues", ev = Eigenvalues[L, 6]];
check["poisson: 6 largest eigenvalues", Round[10^6 Total[ev]]];
bench["poisson: 6 smallest eigenvalues", ev = Eigenvalues[L, -6]];
check["poisson: 6 smallest eigenvalues", Round[10^6 Total[ev]]];

(* ---- sparse eigensystem / singular values with a count (ARPACK) ----------
   Last in the file on purpose: these are the cases still being built, and a
   timeout here costs only these rows.

   NOT the bare Laplacian L.  Its extreme eigenvalues are separated by ~3e-4
   relative, and at that gap the reference's default Arnoldi stops early: the
   two pre-existing L cases above read 47985493 / 4684 there against the exact
   47995643 / 4357 (4 - 2 cos(i pi/301) - 2 cos(j pi/301), which scipy and the
   closed form agree on).  H = L + V adds the potential
   V(i, j) = 0.001 (i^2 + 2 j^2), which separates the ends of the spectrum, so
   all three systems converge to the same six.  Likewise the singular values
   are of a graded bidiagonal-plus-band G (sigma ~ 100/k), not of Dm, whose
   top singular value 91.398 has multiplicity >= 8 (the reference returns five
   copies and a 6th from below it).  The vectors are unit-normalised, so the
   second Eigensystem check element is 1000 x 6. *)

Hpot = L + SparseArray[Band[{1, 1}] -> 0.001 Flatten[Table[i^2 + 2 j^2, {i, $m}, {j, $m}]], {$m^2, $m^2}];
bench["schrodinger: Eigensystem[L + V, 6]", es = Eigensystem[Hpot, 6], 1];
check["schrodinger: Eigensystem[L + V, 6]",
      {Round[10^6 Total[First[es]]], Round[1000 Total[Flatten[Abs[Last[es]]^2]]]}];
$g = 20000;
G = SparseArray[{Band[{1, 1}] -> 100./Range[$g], Band[{1, 2}] -> 0.1, Band[{4, 1}] -> 0.2}, {$g, $g}];
bench["graded 20000: SingularValueList[G, 6]", sv = SingularValueList[G, 6], 1];
check["graded 20000: SingularValueList[G, 6]", Round[10^6 Total[sv]]];
