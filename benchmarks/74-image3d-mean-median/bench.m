v = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 251]]/251, {z, 32}, {y, 48}, {x, 64}], "Real"];
Print["voxels ", 32*48*64];
MeanFilter[v, 2];
Do[Print["MeanFilter   r=", r, "  ", 1000.*First[AbsoluteTiming[Do[MeanFilter[v, r], {5}]]]/5, " ms"], {r, {1, 2, 4}}];
Do[Print["MedianFilter r=", r, "  ", 1000.*First[AbsoluteTiming[Do[MedianFilter[v, r], {3}]]]/3,
   " ms   (window ", (2r+1)^3, ")"], {r, {1, 2}}];
