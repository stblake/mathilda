(* ==========================================================================
   Experiment 53 -- Matrix decompositions (group E, advanced numerical analysis)
   ==========================================================================
   WHAT IT MEASURES.  The dense-decomposition drivers in src/linalg/ --
   ludecomp.c, qrdecomp.c, svdecomp.c, inv.c (PseudoInverse), matrank.c,
   nullspace.c -- all LAPACK-accelerated on the machine-real path.  Group B
   (18-dense-linalg-drivers) times Dot/Solve/Inverse/Det; this goes one level
   down to the factorisations themselves.

   Like group B, EXECUTION not algorithm: numpy/scipy call the same class of
   LAPACK routine, so a gap here is overhead or a missing fast path.

   TIMING runs on large random matrices; CHECKS run on the small exactly
   reproducible spdIntMatrix[8] (from data.m) or a fixed integer matrix, so all
   three systems agree.  Decomposition outputs are non-unique, so every check is
   a convention-invariant scalar: SVD -> sum of singular values, QR -> sum of
   |diag(R)|, LU -> |det| via product of the combined-factor diagonal, rank and
   nullity -> exact integers, PseudoInverse -> total (it is unique).

   CholeskyDecomposition is ABSENT in Mathilda (potrf is internal-only); the
   benchIf case records that as a feature gap while scipy still measures it.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"SingularValueDecomposition", "QRDecomposition", "LUDecomposition",
         "PseudoInverse", "MatrixRank", "NullSpace", "CholeskyDecomposition"}];

seed[];
n  = 400;
a  = rand01[{n, n}];                       (* general square, timing *)
r  = rand01[{300, 150}];                   (* tall, for PseudoInverse timing *)
w  = rand01[{100, 200}];                   (* wide, for NullSpace timing *)
spd = With[{m = rand01[{n, n}]}, Transpose[m].m + n IdentityMatrix[n]];  (* SPD, Cholesky timing *)

sm      = N[spdIntMatrix[8]];              (* check matrix: well-conditioned *)
rankMat = Table[i + j, {i, 6}, {j, 6}];    (* exact rank 2 (outer sum) *)
rect    = N[Table[1/(i + j - 1), {i, 6}, {j, 4}]];  (* full column rank, pinv well-defined *)

bench["SingularValueDecomposition 400x400", SingularValueDecomposition[a];, 1];
check["SingularValueDecomposition 400x400",
      Round[10^6 Total[Diagonal[SingularValueDecomposition[sm][[2]]]]]];

bench["QRDecomposition 400x400", QRDecomposition[a];];
check["QRDecomposition 400x400",
      Round[10^6 Total[Abs[Diagonal[QRDecomposition[sm][[2]]]]]]];

bench["LUDecomposition 400x400", LUDecomposition[a];];
check["LUDecomposition 400x400",
      Round[10^4 Abs[Times @@ Diagonal[LUDecomposition[sm][[1]]]]]];

bench["PseudoInverse 300x150", PseudoInverse[r];];
check["PseudoInverse 300x150", Round[10^6 Total[Flatten[PseudoInverse[rect]]]]];

bench["MatrixRank 400x400", MatrixRank[a];];
check["MatrixRank 400x400", MatrixRank[rankMat]];

(* NullSpace of a float matrix now takes the LAPACK SVD machine path
   (dgesdd, ~6 ms; was seconds through the symbolic row reducer).  Still
   benchOnce -- one run, no warm-up -- so the harness stays comparable. *)
benchOnce["NullSpace 100x200", NullSpace[w];];
check["NullSpace 100x200", Length[NullSpace[rankMat]]];

benchIf["CholeskyDecomposition 400x400", "CholeskyDecomposition",
        CholeskyDecomposition[spd];, 1];
checkIf["CholeskyDecomposition 400x400", "CholeskyDecomposition",
        Round[10^6 Total[Diagonal[CholeskyDecomposition[sm]]]]];
