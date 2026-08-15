v = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 251]]/251, {z, 32}, {y, 48}, {x, 64}], "Real"];
Print["voxels ", 32*48*64];
Dilation[v, 2];
Do[Print["Dilation r=", r, "  ", 1000.*First[AbsoluteTiming[Do[Dilation[v, r], {5}]]]/5,
   " ms   (box is ", (2r+1)^3, " voxels)"], {r, {1, 2, 4, 8}}];
Print["Opening  r=4  ", 1000.*First[AbsoluteTiming[Do[Opening[v, 4], {5}]]]/5, " ms"];
Print["Closing  r=4  ", 1000.*First[AbsoluteTiming[Do[Closing[v, 4], {5}]]]/5, " ms"];
