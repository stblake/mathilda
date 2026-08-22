(* Experiment 90 -- NMinimize robustness testbed (full corpus).
   Companion to 89-nminimize-nmaximize, mirroring the 79/80, 81/82, 83/84 and
   85/86 experiment+testbed pairing.

   89 applies a fair-comparison envelope so every ratio compares like with
   like.  THIS file is where the excluded landscapes go, and it asks a
   different question.

   THE CHECK CARRIES SOLUTION QUALITY, NOT THE OBJECTIVE VALUE.  Each case
   emits Boole[Abs[fbest - fstar] < tol] against the published global optimum:
   1 = found it, 0 = did not.  So a CHECK-FAIL here means exactly "the systems
   disagree about whether they solved the problem", which is the robustness
   gap made legible, rather than a discarded row nobody can interpret.  A row
   where every system emits 1 still gets a legitimate timing.

   This is the only way the harness can express robustness: run_all.py
   discards the timings of any CHECK-FAIL row (run_all.py:521-523), so
   encoding quality INTO the check is what turns a thrown-away row into a
   measurement. *)

Get["../harness.m"];

require[{"NMinimize", "NMaximize"}];

(* solved[f, fstar] -> 1 if the run reached the published optimum. 1e-3 is
   loose on purpose: these are the hard landscapes, and the question is
   "right basin or not", not "how many digits". *)
solved[f_, fstar_] := Boole[Abs[f - fstar] < 0.001];

(* ---- Schwefel 5-D: excluded from 89, Mathilda lands far away ----------- *)
bench["T1 schwefel 5d",
  NMinimize[2094.9144517 - (a Sin[Sqrt[Abs[a]]] + b Sin[Sqrt[Abs[b]]] + c Sin[Sqrt[Abs[c]]] +
    d Sin[Sqrt[Abs[d]]] + e Sin[Sqrt[Abs[e]]]), {a, b, c, d, e}]];
check["T1 schwefel 5d",
  solved[First[NMinimize[2094.9144517 - (a Sin[Sqrt[Abs[a]]] + b Sin[Sqrt[Abs[b]]] + c Sin[Sqrt[Abs[c]]] +
    d Sin[Sqrt[Abs[d]]] + e Sin[Sqrt[Abs[e]]]), {a, b, c, d, e}]], 0]];

(* ---- Griewank 5-D: excluded from 89, Mathilda returns 0.0246 ----------- *)
bench["T2 griewank 5d",
  NMinimize[1 + (a^2+b^2+c^2+d^2+e^2)/4000 -
    Cos[a] Cos[b/Sqrt[2]] Cos[c/Sqrt[3]] Cos[d/2] Cos[e/Sqrt[5]], {a,b,c,d,e}]];
check["T2 griewank 5d",
  solved[First[NMinimize[1 + (a^2+b^2+c^2+d^2+e^2)/4000 -
    Cos[a] Cos[b/Sqrt[2]] Cos[c/Sqrt[3]] Cos[d/2] Cos[e/Sqrt[5]], {a,b,c,d,e}]], 0]];

(* ---- Drop-wave: moved here from 89 because SCIPY is the one that fails -- *)
bench["T3 drop-wave 2d",
  NMinimize[-(1 + Cos[12 Sqrt[x^2 + y^2]])/(0.5 (x^2 + y^2) + 2), {x, y}]];
check["T3 drop-wave 2d",
  solved[First[NMinimize[-(1 + Cos[12 Sqrt[x^2 + y^2]])/(0.5 (x^2 + y^2) + 2), {x, y}]], -1]];

(* ---- Rastrigin 10-D: the dimension where DE starts to strain ----------- *)
bench["T4 rastrigin 10d",
  NMinimize[100 + (a^2-10 Cos[2 Pi a]) + (b^2-10 Cos[2 Pi b]) + (c^2-10 Cos[2 Pi c]) +
    (d^2-10 Cos[2 Pi d]) + (e^2-10 Cos[2 Pi e]) + (f^2-10 Cos[2 Pi f]) +
    (g^2-10 Cos[2 Pi g]) + (h^2-10 Cos[2 Pi h]) + (k^2-10 Cos[2 Pi k]) +
    (m^2-10 Cos[2 Pi m]), {a,b,c,d,e,f,g,h,k,m}]];
check["T4 rastrigin 10d",
  solved[First[NMinimize[100 + (a^2-10 Cos[2 Pi a]) + (b^2-10 Cos[2 Pi b]) + (c^2-10 Cos[2 Pi c]) +
    (d^2-10 Cos[2 Pi d]) + (e^2-10 Cos[2 Pi e]) + (f^2-10 Cos[2 Pi f]) +
    (g^2-10 Cos[2 Pi g]) + (h^2-10 Cos[2 Pi h]) + (k^2-10 Cos[2 Pi k]) +
    (m^2-10 Cos[2 Pi m]), {a,b,c,d,e,f,g,h,k,m}]], 0]];

(* ---- Styblinski-Tang 5-D ----------------------------------------------- *)
bench["T5 styblinski-tang 5d",
  NMinimize[(1/2)((a^4-16a^2+5a)+(b^4-16b^2+5b)+(c^4-16c^2+5c)+(d^4-16d^2+5d)+(e^4-16e^2+5e)),
    {a,b,c,d,e}]];
check["T5 styblinski-tang 5d",
  solved[First[NMinimize[(1/2)((a^4-16a^2+5a)+(b^4-16b^2+5b)+(c^4-16c^2+5c)+(d^4-16d^2+5d)+(e^4-16e^2+5e)),
    {a,b,c,d,e}]], -195.830828518]];

(* ---- Bukin N.6: a razor ridge, notoriously hard ------------------------ *)
bench["T6 bukin n6",
  NMinimize[100 Sqrt[Abs[y - 0.01 x^2]] + 0.01 Abs[x + 10], {x, y}]];
check["T6 bukin n6",
  solved[First[NMinimize[100 Sqrt[Abs[y - 0.01 x^2]] + 0.01 Abs[x + 10], {x, y}]], 0]];

(* ---- Eggholder 2-D ------------------------------------------------------ *)
bench["T7 eggholder 2d",
  NMinimize[-(y + 47) Sin[Sqrt[Abs[y + x/2 + 47]]] - x Sin[Sqrt[Abs[x - (y + 47)]]], {x, y}]];
check["T7 eggholder 2d",
  solved[First[NMinimize[-(y + 47) Sin[Sqrt[Abs[y + x/2 + 47]]] - x Sin[Sqrt[Abs[x - (y + 47)]]], {x, y}]],
    -959.6406627]];

(* ---- F*: the constraint-feasibility gap that 89's C* block rounds away -- *)
(* 89 checks constrained cases at 10^3 so the speed race survives Mathilda's
   ~1e-4 infeasibility.  Here the check IS the feasibility: does the returned
   point actually satisfy x + y >= 2?  Mathilda returns x + y = 1.99990, which
   does not.  Recording it here means the defect is measured rather than
   rounded away and forgotten. *)
F1sol = {x, y} /. Last[NMinimize[{x^2 + y^2, x + y >= 2}, {x, y}]];
bench["F1 ineq feasibility", NMinimize[{x^2 + y^2, x + y >= 2}, {x, y}]];
check["F1 ineq feasibility", Boole[(First[F1sol] + Last[F1sol]) >= 2]];

F2sol = {x, y} /. Last[NMinimize[{x^2 + y^2, x + y == 2}, {x, y}]];
bench["F2 eq feasibility", NMinimize[{x^2 + y^2, x + y == 2}, {x, y}]];
check["F2 eq feasibility", Boole[Abs[First[F2sol] + Last[F2sol] - 2] < 1.*^-9]];
