(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 240 -- SIAM 100-Digit Challenge, sealed holdout B for Appendix D
   problem 5 (min over monic p of ||p(A)||_2, A the 48 x 48 Grcar Toeplitz
   matrix): ten other minimisations of a matrix-norm objective over the
   coefficients of a monic polynomial in a fixed matrix.

     bd5.1   A = diag(Cos[j Pi/6], j = 0..6), deg 6, 2-norm       = 1/32
     bd5.2   A = Grcar 48, deg 8, FROBENIUS norm (smooth)          = 6785.257002...
     bd5.3   A = Grcar 48, deg 4, 2-norm                            = 45.55617241...
     bd5.4   A = 11 x 11 Kac (Clement) matrix, deg 4, 2-norm       = 1152
     bd5.5   A = tridiag(-1, 2, -1), n = 20, deg 4, 2-norm          = 1.873659770...
     bd5.6   A = tridiagonal Toeplitz (sub 1, super 1/4), n = 40, deg 3, 2-norm
                                                                    = 1.031656614...
     bd5.7   A = Grcar 48, deg 1: min_c ||A + c I||_2               = 3.129559188...
     bd5.8   A = Grcar 24, deg 6, 2-norm                            = 272.0889162...
     bd5.9   A = HilbertMatrix[8], deg 3, 2-norm                    = 0.005932729334...
     bd5.10  A = Toeplitz n = 32, 1 on the first subdiagonal and the first two
             superdiagonals, 0 on the diagonal; deg 3, 2-norm      = 7.231703771...

   Grcar n: -1 on the first subdiagonal, 1 on the diagonal and the first three
   superdiagonals.  Kac n = 11: zero diagonal, off-diagonal Sqrt[j (11 - j)]
   (eigenvalues -10, -8, ..., 10).  For a SYMMETRIC A (bd5.1, 4, 5, 9) the
   2-norm objective is the discrete Chebyshev problem max_i |p(lambda_i)|,
   whose optimum is a crease where k+1 eigenvalues of p(A) are equal in
   modulus (equioscillation).  bd5.8 is the nonnormal analogue of experiment
   180: at the optimum the top singular value is double, the objective is
   not differentiable there, and a smooth optimiser stalls on the crease
   (scipy BFGS 272.167, Nelder-Mead 463).  bd5.3/6/7/10 are nonnormal too,
   but their optimum happens to be a SMOOTH point (simple top singular value)
   -- a correct generic minimiser gets them; bd5.2 is the Frobenius norm, a
   least-squares problem in disguise.  The easy rows are part of the probe:
   a system that only handles the smooth ones has not fixed the crease.

   Idiomatic statement, as in experiment 180: FindMinimum of Norm[p(A), .]
   over the k lower coefficients from c = 0 (p = x^k), the powers A^j
   precomputed exactly (glue); machine precision and WorkingPrecision -> 60.
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindMinimum", "Norm", "MatrixPower", "HilbertMatrix"}];

toep[n_, rules_] := Table[Replace[j - i, Append[rules, _ -> 0]], {i, n}, {j, n}];
grcar[n_] := toep[n, {-1 -> -1, 0 -> 1, 1 -> 1, 2 -> 1, 3 -> 1}];
kac = Table[If[Abs[i - j] == 1, Sqrt[Min[i, j] (11 - Min[i, j])], 0], {i, 11}, {j, 11}];

(* one case: monic degree k in A, norm nrm; start c = 0 *)
mcase[lab_String, A_, k_Integer, nrm_] :=
  (* A^0 is IdentityMatrix: wolframscript's MatrixPower[A, 0] refuses a
     singular A (MatrixPower::sing -- bd5.1's diagonal has a 0, Kac's 11x11
     is singular) *)
  Module[{pw = Prepend[Table[MatrixPower[A, j], {j, 1, k}], IdentityMatrix[Length[A]]], f, cv},
    f[c__?NumericQ] := Norm[pw[[k + 1]] + Total[{c} pw[[1 ;; k]]], nrm];
    cv = Table[cc[j], {j, k}];
    vrun[lab <> " machine", First[FindMinimum[f @@ cv, Transpose[{cv, ConstantArray[0., k]}]]],
         10, 1, 60];
    vrun[lab <> " wp50", First[FindMinimum[f @@ cv, Transpose[{cv, ConstantArray[0, k]}],
                                         WorkingPrecision -> 60]], 48, 1, 120]];

mcase["bd5.1", DiagonalMatrix[Table[Cos[j Pi/6], {j, 0, 6}]], 6, 2];
mcase["bd5.2", grcar[48], 8, "Frobenius"];
mcase["bd5.3", grcar[48], 4, 2];
mcase["bd5.4", kac, 4, 2];
mcase["bd5.5", toep[20, {-1 -> -1, 0 -> 2, 1 -> -1}], 4, 2];
mcase["bd5.6", toep[40, {-1 -> 1, 1 -> 1/4}], 3, 2];
mcase["bd5.7", grcar[48], 1, 2];
mcase["bd5.8", grcar[24], 6, 2];
mcase["bd5.9", HilbertMatrix[8], 3, 2];
mcase["bd5.10", toep[32, {-1 -> 1, 1 -> 1, 2 -> 1}], 3, 2];
