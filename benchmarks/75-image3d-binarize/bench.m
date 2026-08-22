v = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 251]]/251, {z, 32}, {y, 48}, {x, 64}], "Real"];
Print["voxels ", 32*48*64];
Binarize[v];
Print["Binarize (Otsu)      ", 1000.*First[AbsoluteTiming[Do[Binarize[v], {5}]]]/5, " ms"];
Do[Print["LocalAdaptive r=", r, "     ", 1000.*First[AbsoluteTiming[Do[LocalAdaptiveBinarize[v, r], {5}]]]/5, " ms"], {r, {2, 4, 8}}];
Print["LocalAdaptive sauvola ", 1000.*First[AbsoluteTiming[Do[LocalAdaptiveBinarize[v, 4, {1, -0.2, 0.}], {5}]]]/5, " ms"];
