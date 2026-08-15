img = Image[Table[N[Mod[i*7 + j*13, 251]]/251, {i, 512}, {j, 512}], "Real"];
t32 = ImageData[img][[100 ;; 131, 200 ;; 231]];
t8  = ImageData[img][[100 ;; 107, 200 ;; 207]];
Do[ImageCorrelate[img, t32, "NormalizedCrossCorrelation"], {2}];
Print["NCC 32x32 template ", 1000.*First[AbsoluteTiming[
   Do[ImageCorrelate[img, t32, "NormalizedCrossCorrelation"], {5}]]]/5, " ms"];
Print["NCC 8x8 template   ", 1000.*First[AbsoluteTiming[
   Do[ImageCorrelate[img, t8, "NormalizedCrossCorrelation"], {5}]]]/5, " ms"];
Print["plain correlate 32 ", 1000.*First[AbsoluteTiming[
   Do[ImageCorrelate[img, t32], {5}]]]/5, " ms"];
