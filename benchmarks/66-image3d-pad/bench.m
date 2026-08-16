v = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 251]]/251, {z, 64}, {y, 96}, {x, 128}], "Real"];
Do[ImagePad[v, 8], {2}];
Print["volume         ", ImageDimensions[v], "  = ", 64*96*128, " voxels"];
Print["ImagePad 8     ", 1000.*First[AbsoluteTiming[Do[ImagePad[v, 8], {10}]]]/10, " ms"];
Print["ImagePad refl  ", 1000.*First[AbsoluteTiming[Do[ImagePad[v, 8, "Reflected"], {10}]]]/10, " ms"];
big = ImagePad[v, 8];
Print["ImageCrop      ", 1000.*First[AbsoluteTiming[Do[ImageCrop[big, {128, 96, 64}], {10}]]]/10, " ms"];
