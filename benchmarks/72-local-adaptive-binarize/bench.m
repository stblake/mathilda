img = Image[Table[N[Mod[i*7 + j*13, 251]]/251, {i, 512}, {j, 512}], "Real"];
LocalAdaptiveBinarize[img, 8];
Do[Print["r=", r, "  ", 1000.*First[AbsoluteTiming[Do[LocalAdaptiveBinarize[img, r], {5}]]]/5, " ms"],
 {r, {2, 8, 16, 32}}];
Print["Binarize (global Otsu)  ", 1000.*First[AbsoluteTiming[Do[Binarize[img], {5}]]]/5, " ms"];
