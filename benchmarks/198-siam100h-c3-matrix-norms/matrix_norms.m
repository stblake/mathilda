(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 198 -- SIAM 100-Digit Challenge, holdout A for Problem 3
   (family hc3: the 2-norm, i.e. the largest singular value, of a matrix
   given by a formula for its entries).  The statement is C3's (experiment
   168): Norm of the n x m section, built directly in the working precision,
   in doubles and at 60 digits.  Infinite operators are chosen so that the
   stated section carries their norm to > 50 digits; the others are the
   finite matrices themselves.

     k   a_ij                               section    scored against
     1   2^-i 3^-j                          100 x 100  ||A|| = 1/sqrt(24)      (rank one, easy)
     2   1/(i! j!)                          30 x 30    ||A|| = I_0(2) - 1      (rank one, easy)
     3   1/(i+j-1)!                         40 x 40    ||A|| (infinite operator)
     4   1/((i+j-1) 2^(i+j))                100 x 100  ||A|| (infinite operator)
     5   1/(i+j-1)        (Hilbert)         100 x 100  ||H_100||
     6   min(i,j)/max(i,j) (Lehmer)         100 x 100  ||L_100||
     7   1/((i+j-1)(i+j-2)/2 + j)^2         100 x 100  ||.|| (C3's entries squared)
     8   1/(i + j^2)                        120 x 60   ||.|| (rectangular)
     9   1/(1+|i-j|)^2    (Toeplitz)        150 x 150  ||.|| (sigma_2/sigma_1 = 0.974: clustered)
     10  cos(i j)/(i+j)                     100 x 100  ||.|| (sign-scrambled, sigma_2/sigma_1 = 0.78)

   True values: truth.tsv; references siam100_refs/hc3_ref.{py,wl}.
   Maxima's SVD is LAPACK dgesvd (machine only; wp50 SKIP); FriCAS has no
   numerical SVD (SKIP).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"Norm"}];

ent[1] = Function[{i, j}, 2^-i 3^-j];
ent[2] = Function[{i, j}, 1/(i! j!)];
ent[3] = Function[{i, j}, 1/(i + j - 1)!];
ent[4] = Function[{i, j}, 1/((i + j - 1) 2^(i + j))];
ent[5] = Function[{i, j}, 1/(i + j - 1)];
ent[6] = Function[{i, j}, Min[i, j]/Max[i, j]];
ent[7] = Function[{i, j}, 1/((i + j - 1) (i + j - 2)/2 + j)^2];
ent[8] = Function[{i, j}, 1/(i + j^2)];
ent[9] = Function[{i, j}, 1/(1 + Abs[i - j])^2];
ent[10] = Function[{i, j}, Cos[i j]/(i + j)];
size = {{100, 100}, {30, 30}, {40, 40}, {100, 100}, {100, 100}, {100, 100},
        {100, 100}, {120, 60}, {150, 150}, {100, 100}};

sec[k_, prec_] := Table[N[ent[k][i, j], prec], {i, size[[k, 1]]}, {j, size[[k, 2]]}];

Do[
  vrun["hc3." <> ToString[k] <> " machine", Norm[sec[k, MachinePrecision]], 10, 3, 60];
  vrun["hc3." <> ToString[k] <> " wp50", Norm[sec[k, 60]], 48, 1, 90],
  {k, 10}];
