img = Image[Table[N[Mod[i*7 + j*13, 251]]/251, {i, 512}, {j, 512}], "Real"];
CornerFilter[img];
Print["CornerFilter MinimumEigenvalue r=2 ", 1000.*First[AbsoluteTiming[Do[CornerFilter[img], {5}]]]/5, " ms"];
Print["CornerFilter Harris r=2            ", 1000.*First[AbsoluteTiming[Do[CornerFilter[img, 2, "Harris"], {5}]]]/5, " ms"];
Print["ImageCorners r=2                   ", 1000.*First[AbsoluteTiming[Do[ImageCorners[img], {3}]]]/3, " ms"];
