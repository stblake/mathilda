(* c3 reference, part 1: sigma_max(A_n) of the n x n section to 120 digits by
   power iteration on A^T A at 160-digit working precision (sigma2/sigma1 ~ 0.11,
   so ~1.9 digits per sweep), for a ladder of n.  Output: c3_trunc.m *)
a[i_, j_] := 1/((i + j - 1) (i + j - 2)/2 + j);
smax[n_] := Module[{A = N[Table[a[i, j], {i, n}, {j, n}], 160], At, v, w, lam = 0, lam0 = -1, k = 0},
  At = Transpose[A]; v = ConstantArray[N[1, 160], n];
  While[k < 400 && Abs[lam - lam0] > 10^-140, lam0 = lam; w = At.(A.v);
    lam = Sqrt[w.w]/Sqrt[v.v]; v = w/Sqrt[w.w]; k++];
  w = At.(A.v); Sqrt[v.w/v.v]];
ns = Join[Range[40, 400, 20], {450, 500, 600, 700, 800}];
res = Table[{n, smax[n]}, {n, ns}];
Put[res, "c3_trunc.m"];
Print[N[res[[-1]], 30]];
