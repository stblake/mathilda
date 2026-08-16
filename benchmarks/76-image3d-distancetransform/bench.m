v = Image3D[Table[If[Mod[z*5 + y*3 + x*7, 23] == 0, 0., 1.], {z, 32}, {y, 48}, {x, 64}], "Real"];
Print["voxels ", 32*48*64, "  background fraction ~", N[1/23]];
DistanceTransform[v];
Print["DistanceTransform  ", 1000.*First[AbsoluteTiming[Do[DistanceTransform[v], {5}]]]/5, " ms"];
sparse = Image3D[Table[If[z == 16 && y == 24 && x == 32, 0., 1.], {z, 32}, {y, 48}, {x, 64}], "Real"];
Print["single seed (worst) ", 1000.*First[AbsoluteTiming[Do[DistanceTransform[sparse], {5}]]]/5, " ms"];
