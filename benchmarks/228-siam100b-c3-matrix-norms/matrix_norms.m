(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 228 -- SIAM 100-Digit Challenge, sealed holdout B for Problem 3
   ("How Far Away Is Infinity?"): the 2-norm (largest singular value) of ten
   more dense matrices -- sections of infinite matrices, some bounded
   operators whose sections converge superexponentially, some not.

     bc3.1   Hilbert a_ij = 1/(i+j-1)                     m: n = 500   w: n = 60
     bc3.2   Hankel a_ij = 1/(i+j-1)!   (||A|| itself:     both n = 40; the
             the section converges to > 70 digits)
     bc3.3   Toeplitz a_ij = 1/(1+(i-j)^2) (Lorentzian)     m: 400    w: 60
     bc3.4   Cesaro a_ij = 1/i (j <= i), 0 (j > i)          m: 500    w: 80
     bc3.5   Toeplitz a_ij = 2^-(i-j) (i >= j), (-1/3)^(j-i) (j > i)
                                                            m: 400    w: 60
     bc3.6   the Problem-3 matrix, RECTANGULAR sections     m: 500x150 w: 100x30
     bc3.7   complex a_jk = exp(2 pi i jk/7)/(j+k)          m: 300    w: 50
     bc3.8   Frank matrix f_ij = n+1-max(i,j) (j >= i-1), 0 below
                                                            m: 300    w: 60
     bc3.9   a_ij = 1/(2^i + 3^j)  (||A|| itself: both n = 150; converged)
     bc3.10  a_ij = 1/sqrt(i+j-1) (unbounded; sections)    m: 1000   w: 120

   Idiomatic statement: Norm of the section, built directly in the working
   precision (machine doubles; 60 digits).  Each row is scored against the
   largest singular value of ITS section (bc3.2, bc3.9: equal to ||A|| to
   > 70 digits).  True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"Norm"}];

mat[1, n_, p_] := Table[N[1/(i + j - 1), p], {i, n}, {j, n}];
mat[2, n_, p_] := Table[N[1/(i + j - 1)!, p], {i, n}, {j, n}];
mat[3, n_, p_] := Table[N[1/(1 + (i - j)^2), p], {i, n}, {j, n}];
mat[4, n_, p_] := Table[N[If[j <= i, 1/i, 0], p], {i, n}, {j, n}];
mat[5, n_, p_] := Table[N[If[i >= j, 2^-(i - j), (-1/3)^(j - i)], p], {i, n}, {j, n}];
mat[6, {m_, n_}, p_] := Table[N[1/((i + j - 1) (i + j - 2)/2 + j), p], {i, m}, {j, n}];
mat[7, n_, p_] := Table[N[Exp[2 Pi I Mod[j k, 7]/7]/(j + k), p], {j, n}, {k, n}];
mat[8, n_, p_] := Table[N[If[j >= i - 1, n + 1 - Max[i, j], 0], p], {i, n}, {j, n}];
mat[9, n_, p_] := Table[N[1/(2^i + 3^j), p], {i, n}, {j, n}];
mat[10, n_, p_] := Table[N[1/Sqrt[i + j - 1], p], {i, n}, {j, n}];

sizes = {{500, 60}, {40, 40}, {400, 60}, {500, 80}, {400, 60},
         {{500, 150}, {100, 30}}, {300, 50}, {300, 60}, {150, 150}, {1000, 120}};

Do[With[{k = k, nm = sizes[[k, 1]], nw = sizes[[k, 2]]},
    vrun["bc3." <> ToString[k] <> " machine", Norm[mat[k, nm, MachinePrecision]], 10, 3, 60];
    vrun["bc3." <> ToString[k] <> " wp50", Norm[mat[k, nw, 60]], 48, 1, 120]],
  {k, 10}];
