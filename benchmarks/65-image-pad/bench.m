img = Image[Table[N[Mod[i*7 + j*13, 251]]/251, {i, 512}, {j, 512}], "Real"];
Do[ImagePad[img, 16], {3}];
t1 = First[AbsoluteTiming[Do[ImagePad[img, 16], {20}]]]/20;
Print["ImagePad 16 (const)  ", 1000.*t1, " ms"];
t2 = First[AbsoluteTiming[Do[ImagePad[img, 16, "Reflected"], {20}]]]/20;
Print["ImagePad 16 reflect  ", 1000.*t2, " ms"];
big = ImagePad[img, 16];
t3 = First[AbsoluteTiming[Do[ImageCrop[big, {512, 512}], {20}]]]/20;
Print["ImageCrop centred    ", 1000.*t3, " ms"];
