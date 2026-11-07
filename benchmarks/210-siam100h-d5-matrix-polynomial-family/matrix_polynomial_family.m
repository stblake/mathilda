(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 210 -- SIAM 100-digit HOLDOUT A, family d5 (best polynomial of
   a matrix; original: experiment 180, min over monic p of degree 8 of
   ||p(A)||_2 for the 48 x 48 Grcar-type Toeplitz matrix).

   min over normalised p of degree k of ||p(A)||:
     hd5.1   A = diag(Cos[j Pi/8], j = 0..8),  monic, k = 8,  2-norm  (= 2^-7)
     hd5.2   A = diag(1, ..., 10),             monic, k = 4,  2-norm
     hd5.3   A = the 48 x 48 matrix of exp. 180, monic, k = 8, FROBENIUS norm
     hd5.4   A = tridiag(-1, 2, -1), n = 20,    monic, k = 5,  2-norm
     hd5.5   A = Grcar(n = 32, 3 superdiagonals) [-1 sub, 1 diag, 1 super x3],
                                               monic, k = 6,  2-norm
     hd5.6   A = Grcar(n = 40, 2 superdiagonals), monic, k = 6, 2-norm
     hd5.7   A = diag(1, ..., 10),  p(0) = 1 (ideal GMRES),  k = 3,  2-norm
     hd5.8   A = upper bidiagonal Toeplitz, 1 on the diagonal, 2 above, n = 12,
                                               monic, k = 4,  2-norm
     hd5.9   A = tridiag(1/4, 0, 1) (nonsymmetric Toeplitz), n = 24,
                                               monic, k = 6,  2-norm
     hd5.10  A = Grcar(n = 64, 3 superdiagonals), monic, k = 10, 2-norm (harder)
   Generalised along the original's axis -- a convex but NONSMOOTH min-max
   (the 2-norm is a largest singular value, which is multiple at the
   minimiser).  Varied in kind: normal matrices, where the answer is a
   discrete Chebyshev (minimax) problem on the spectrum and the crease has
   k+1 active points (hd5.1, .2, .4, .7); non-normal Toeplitz matrices with a
   two-fold crease (hd5.5, .6, .8, .9, .10); the smooth Frobenius-norm
   relative (hd5.3, a least-squares problem); and the p(0) = 1 normalisation
   of ideal GMRES (hd5.7).

   Idiomatic statement, as experiment 180: FindMinimum of Norm[p(A)] (or
   Norm[p(A), "Frobenius"]) over the k free coefficients from c = 0, the
   powers A^j precomputed exactly (glue), machine precision and
   WorkingPrecision -> 60.  References: siam100_refs/hd5_ref.py and .wl (exact least
   squares / alternation systems / KKT Newton on the crease); truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindMinimum", "Norm", "MatrixPower"}];

toep[n_, d_] := Table[Lookup[d, j - i, 0], {i, n}, {j, n}];
grcar[n_, s_] := toep[n, Join[<|-1 -> -1, 0 -> 1|>, AssociationThread[Range[s] -> 1]]];
(* {A, k, normalisation, norm} *)
probs = {
  {DiagonalMatrix[Table[Cos[j Pi/8], {j, 0, 8}]], 8, "monic", 2},
  {DiagonalMatrix[Range[10]], 4, "monic", 2},
  {grcar[48, 3], 8, "monic", "Frobenius"},
  {toep[20, <|-1 -> -1, 0 -> 2, 1 -> -1|>], 5, "monic", 2},
  {grcar[32, 3], 6, "monic", 2},
  {grcar[40, 2], 6, "monic", 2},
  {DiagonalMatrix[Range[10]], 3, "gmres", 2},
  {toep[12, <|0 -> 1, 1 -> 2|>], 4, "monic", 2},
  {toep[24, <|-1 -> 1/4, 1 -> 1|>], 6, "monic", 2},
  {grcar[64, 3], 10, "monic", 2}};

(* B0 + Sum_j c_j B_j : monic -> B0 = A^k, B_j = A^(j-1); gmres -> B0 = I, B_j = A^j *)
Do[With[{A = probs[[i, 1]], k = probs[[i, 2]], nz = probs[[i, 3]], nm = probs[[i, 4]]},
   Module[{pw = Table[MatrixPower[A, j], {j, 0, k}], b0, bs, obj, obj60, vars, ex},
     If[nz === "monic", b0 = pw[[k + 1]]; bs = pw[[1 ;; k]], b0 = pw[[1]]; bs = pw[[2 ;; k + 1]]];
     (* exact rational powers as in exp. 180; the algebraic Cos[j Pi/8] powers
        of hd5.1 are numericised once (glue), at each tier's precision *)
     ex = MatrixQ[A, Head[#] === Integer || Head[#] === Rational &];
     vars = Table[cc[i, j], {j, k}];
     With[{b0m = If[ex, b0, N[b0]], bsm = If[ex, bs, N[bs]],
           b060 = If[ex, b0, N[b0, 80]], bs60 = If[ex, bs, N[bs, 80]]},
       obj[x__?NumericQ] := Norm[b0m + Total[{x} bsm], nm];
       obj60[x__?NumericQ] := Norm[b060 + Total[{x} bs60], nm]];
     vrun["hd5." <> ToString[i] <> " machine",
       First[FindMinimum[obj @@ vars, Transpose[{vars, ConstantArray[0., k]}]]], 10, 1, 60];
     vrun["hd5." <> ToString[i] <> " wp50",
       First[FindMinimum[obj60 @@ vars, Transpose[{vars, ConstantArray[0, k]}],
         WorkingPrecision -> 60]], 48, 1, If[i == 10, 120, 60]]]],
  {i, 1, 10}];
