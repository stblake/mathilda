(* hc3 (holdout A, family C3) reference, method B: WL SingularValueList
   (largest singular value) of each matrix at precision 120; k = 3, 4 at the
   larger sections n = 60 / 160 (the infinite operator).  Method A:
   hc3_ref.py (mpmath eigsy of A^T A at dps 90; closed forms for k = 1, 2). *)
mat[1] = Table[2^-i 3^-j, {i, 100}, {j, 100}];
mat[2] = Table[1/(i! j!), {i, 30}, {j, 30}];
mat[3] = Table[1/(i + j - 1)!, {i, 60}, {j, 60}];
mat[4] = Table[1/((i + j - 1) 2^(i + j)), {i, 160}, {j, 160}];
mat[5] = Table[1/(i + j - 1), {i, 100}, {j, 100}];
mat[6] = Table[Min[i, j]/Max[i, j], {i, 100}, {j, 100}];
mat[7] = Table[1/((i + j - 1) (i + j - 2)/2 + j)^2, {i, 100}, {j, 100}];
mat[8] = Table[1/(i + j^2), {i, 120}, {j, 60}];
mat[9] = Table[1/(1 + Abs[i - j])^2, {i, 150}, {j, 150}];
mat[10] = Table[Cos[i j]/(i + j), {i, 100}, {j, 100}];
Do[Print[k, " ", ToString[NumberForm[First[SingularValueList[N[mat[k], 120], 1]], 70]]], {k, 10}];
Print["closed 1: ", ToString[NumberForm[N[Sqrt[1/24], 70], 70]]];
Print["closed 2: ", ToString[NumberForm[N[BesselI[0, 2] - 1, 70], 70]]];
