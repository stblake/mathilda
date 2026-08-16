img = Image[Table[N[Mod[i*7 + j*13, 251]]/251, {i, 512}, {j, 512}], "Real"];
(* NON-separable kernels, so the separable path never applies *)
mk[n_] := Table[N[Mod[i*5 + j*3, 11]] - 5, {i, 1, n}, {j, 1, n}];
Do[Module[{k = mk[n], t},
   ImageConvolve[img, k];
   t = First[AbsoluteTiming[Do[ImageConvolve[img, k], {3}]]]/3;
   Print["k=", n, "x", n, "  rank=", MatrixRank[k], "  ", 1000.*t, " ms"]],
 {n, {3, 5, 7, 9, 11, 15, 21, 32}}];
