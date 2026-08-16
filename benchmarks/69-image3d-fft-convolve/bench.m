v = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 251]]/251, {z, 64}, {y, 96}, {x, 128}], "Real"];
mk[n_] := Table[N[Mod[m*5 + i*3 + j*2, 11]] - 5, {m, 1, n}, {i, 1, n}, {j, 1, n}];
Do[Module[{k = mk[n], t},
   ImageConvolve[v, k];
   t = First[AbsoluteTiming[Do[ImageConvolve[v, k], {3}]]]/3;
   Print["k=", n, "^3 (", n^3, " taps)  ", 1000.*t, " ms"]],
 {n, {3, 5, 7, 9}}];
