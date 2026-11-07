(* D1 reference, second implementation (wolframscript; different truncation
   from d1_ref.py; moments and tail solve in 150-digit
   arithmetic -- an exact rational solve of the 90x90 system is too slow):
   Schmelzer-Baillie power-sum recursion for sum 1/n over n avoiding "42".
   Explicit sum over n < 10^K0 (K0 = 3), moments j = 1..J (J = 45), the tail
   e.(I - M)^-1 M psi(K0) by LinearSolve at 150 digits. *)
k0 = 3; jj = 45;
adm = Select[Range[10^k0 - 1], StringFreeQ[IntegerString[#], "42"] &];
head = N[Total[1/adm], 150];
st[n_] := If[Mod[n, 10] == 4, 2, 1];
last = Select[adm, IntegerLength[#] == k0 &];
psi = Flatten[Table[Total[N[Select[last, st[#] == s &], 150]^-j], {s, 1, 2}, {j, 1, jj}]];
idx[s_, j_] := (s - 1) jj + j;
trans = Join[Table[{1, d, If[d == 4, 2, 1]}, {d, 0, 9}],
             Table[{2, d, If[d == 4, 2, 1]}, {d, DeleteCases[Range[0, 9], 2]}]];
m = ConstantArray[0, {2 jj, 2 jj}];
Do[{i, d, l} = tr;
   Do[m[[idx[l, j], idx[i, j + w]]] += (-1)^w Binomial[j + w - 1, w] If[w == 0, 1, d^w]/10^(j + w),
      {j, 1, jj}, {w, 0, jj - j}], {tr, trans}];
y = LinearSolve[N[IdentityMatrix[2 jj] - m, 150], N[m . psi, 150]];
Print[N[head + y[[idx[1, 1]]] + y[[idx[2, 1]]], 75]];
