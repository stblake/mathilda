v = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 251]]/251, {z, 32}, {y, 48}, {x, 64}], "Real"];
ImageReflect[v];
Do[Print["ImageReflect ", sd, "  ", 1000.*First[AbsoluteTiming[Do[ImageReflect[v, sd], {10}]]]/10, " ms"],
 {sd, {Top, Left, Front}}];
