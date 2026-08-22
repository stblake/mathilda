v = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 251]]/251, {z, 32}, {y, 48}, {x, 64}], "Real"];
Print["voxels ", 32*48*64];
GradientFilter[v];
Print["GradientFilter      ", 1000.*First[AbsoluteTiming[Do[GradientFilter[v], {5}]]]/5, " ms"];
Print["DerivativeFilter x  ", 1000.*First[AbsoluteTiming[Do[DerivativeFilter[v, {0, 0, 1}], {5}]]]/5, " ms"];
Print["DerivativeFilter xx ", 1000.*First[AbsoluteTiming[Do[DerivativeFilter[v, {0, 0, 2}], {5}]]]/5, " ms"];
