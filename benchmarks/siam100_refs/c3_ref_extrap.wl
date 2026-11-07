(* c3 reference, part 2: ||A|| = lim sigma_max(A_n).  Polynomial (Richardson)
   extrapolation in h = 1/n of the 120-digit sections from c3_ref_trunc.wl:
   interpolate the last m points by a polynomial in h and evaluate at h = 0;
   the agreement across m and across point subsets is the digit count. *)
res = Get["c3_trunc.m"];
ext[pts_] := InterpolatingPolynomial[{1/#[[1]], #[[2]]} & /@ pts, h] /. h -> 0;
Do[Print[m, "  ", N[ext[Take[res, -m]], 60]], {m, 8, Length[res], 2}];
Print["even-indexed: ", N[ext[res[[2 ;; ;; 2]]], 60]];
Print["odd-indexed:  ", N[ext[res[[1 ;; ;; 2]]], 60]];
