v = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 251]]/251, {z, 32}, {y, 48}, {x, 64}], "Real"];
Print["voxels ", 32*48*64];
CornerFilter[v];
Print["CornerFilter 3D MinimumEigenvalue ", 1000.*First[AbsoluteTiming[Do[CornerFilter[v], {3}]]]/3, " ms"];
Print["CornerFilter 3D Harris            ", 1000.*First[AbsoluteTiming[Do[CornerFilter[v, 2, "Harris"], {3}]]]/3, " ms"];
Print["ImageCorners 3D                   ", 1000.*First[AbsoluteTiming[Do[ImageCorners[v], {3}]]]/3, " ms  (", Length[ImageCorners[v]], " corners)"];
