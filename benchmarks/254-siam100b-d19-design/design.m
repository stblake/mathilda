(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 254 -- SIAM 100-Digit Challenge, SEALED HOLDOUT B, family D19
   (Appendix D Problem 19, optimal design: inf of Int u over conductivities
   c in {1, 100}; exp 194 asks its relaxed finite proxy).

   Ten new finite convex "variable thickness" design problems with the
   capability exp 194 tests: a constrained minimisation (box + one linear
   equality) of a black-box objective built with LinearSolve,
       J(c) = F . U,  K(c) U = F,  K(c) = Sum_e c_e K_e,
       lo <= c_e <= hi,  Sum_e c_e = V.
   J is convex in c, so the minimum VALUE is unique (truth.tsv:
   siam100_refs/bd19_ref.py KKT Newton at 90 digits == bd19_ref.gp exact
   rational duality certificate).
     bd19.1   1-D bar [0,1], 8 P1 elements, u(0) = 0, free end, f = 1,  c in [1,10],  Sum 32   (easy)
     bd19.2   1-D bar, 12 elements, u(0) = u(1) = 0, f = 1,             c in [10,100], Sum 606
     bd19.3   unit square, 4x4 Q1 cells, u = 0 on the boundary, f = 1,  c in [1,10],  Sum 48
     bd19.4   rectangle [0,2]x[0,1], 8x4 Q1 cells, f = 1,               c in [1,100], Sum 1616
     bd19.5   unit square, 6x6 Q1, unit POINT load at the centre node,  c in [1,100], Sum 1818
     bd19.6   unit square, 6x6 Q1, u = 0 ONLY on {y = 0, x <= 1/3} (a
              heat sink; zero flux elsewhere), f = 1,                    c in [1,100], Sum 360
     bd19.7   resistor network: 4x4-node grid graph, 24 edge conductances,
              node (0,0) grounded, unit current into (1,2); J = effective
              resistance,                                                 c in [1,10],  Sum 60
     bd19.8   unit square, 4x4 cells each cut along its (a,b)-(a+1,b+1)
              diagonal into two P1 triangles sharing one c (16 unknowns),
              f = 1,                                                      c in [1,100], Sum 808
     bd19.9   L-shape: 6x6 Q1 grid on the unit square minus the cells
              a, b >= 3 (27 cells), u = 0 on its boundary, f = 1,       c in [1,100], Sum 1350
     bd19.10  unit square, 10x10 Q1 (100 unknowns), f = 1,              c in [1,100], Sum 5050 (harder)
   f = 1: the consistent load of a unit source (each element adds the
   integral of each of its basis functions).  Q1 element stiffness
   c_e (1/6)[[4,-1,-2,-1],[-1,4,-1,-2],[-2,-1,4,-1],[-1,-2,-1,4]] (corners
   counter-clockwise, square cells); P1 bar c_e/h [[1,-1],[-1,1]]; P1
   triangle pair c_e (1/2)[[2,-1,0,-1],[-1,2,-1,0],[0,-1,2,-1],[-1,0,-1,2]];
   resistor c_e [[1,-1],[-1,1]].
   Idiomatic statement (as exp 194): FindMinimum with the bounds and the
   equality, the objective a numeric black box assembled from the element
   list (glue) and solved by LinearSolve.  Rows: "bd19.k machine" and
   "bd19.k wp50" (WorkingPrecision -> 60, PrecisionGoal -> 50).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindMinimum", "LinearSolve"}];

(* A problem is {nDof, {{dofs, Kloc}, ...}, F}; a dof 0 is a fixed (u = 0) node. *)
kq = {{4, -1, -2, -1}, {-1, 4, -1, -2}, {-2, -1, 4, -1}, {-1, -2, -1, 4}}/6;

loadOf[nDof_, els_, wts_] := Module[{F = ConstantArray[0, nDof]},
  Do[Do[If[els[[e, 1, r]] > 0, F[[els[[e, 1, r]]]] += wts[[r]]], {r, Length[wts]}], {e, Length[els]}];
  F];

bar[n_, clampRight_] := Module[{h = 1/n, dof, els},
  dof[i_] := If[i == 0 || (clampRight && i == n), 0, i];
  els = Table[{{dof[e], dof[e + 1]}, {{1, -1}, {-1, 1}}/h}, {e, 0, n - 1}];
  {If[clampRight, n - 1, n], els, loadOf[If[clampRight, n - 1, n], els, {h/2, h/2}]}];

(* Q1 on an nx x ny grid of square cells of side h; cellOk[a, b]; nodeFree[i, j];
   point: Automatic (unit source) or a node {i, j} carrying a unit point load *)
q1[nx_, ny_, h_, cellOk_, nodeFree_, point_ : Automatic] := Module[{cells, used, free, dof, els, F},
  cells = Select[Tuples[{Range[0, nx - 1], Range[0, ny - 1]}], cellOk @@ # &];
  used = Union @@ ({#, # + {1, 0}, # + {1, 1}, # + {0, 1}} & /@ cells);
  free = Select[used, nodeFree @@ # &];
  dof = AssociationThread[free -> Range[Length[free]]];
  els = {Lookup[dof, Key /@ {#, # + {1, 0}, # + {1, 1}, # + {0, 1}}, 0], kq} & /@ cells;
  F = If[point === Automatic, loadOf[Length[free], els, ConstantArray[h^2/4, 4]],
    ReplacePart[ConstantArray[0, Length[free]], dof[point] -> 1]];
  {Length[free], els, F}];

interiorQ[nx_, ny_] := (0 < #1 < nx && 0 < #2 < ny) &;

(* each cell = two P1 right triangles (diagonal (a,b)-(a+1,b+1)) sharing ONE c: the summed
   cell matrix is the 5-point stencil, the load h^2 {1/3, 1/6, 1/3, 1/6} on the corners *)
p1cell[n_] := Module[{h = 1/n, dof, kc, els},
  dof[{i_, j_}] := If[0 < i < n && 0 < j < n, (i - 1) (n - 1) + j, 0];
  kc = {{2, -1, 0, -1}, {-1, 2, -1, 0}, {0, -1, 2, -1}, {-1, 0, -1, 2}}/2;
  els = Flatten[Table[{dof /@ {{a, b}, {a + 1, b}, {a + 1, b + 1}, {a, b + 1}}, kc},
      {a, 0, n - 1}, {b, 0, n - 1}], 1];
  {(n - 1)^2, els, loadOf[(n - 1)^2, els, h^2 {1/3, 1/6, 1/3, 1/6}]}];

gridGraph[m_] := Module[{nodes, dof, els, F},
  nodes = Flatten[Table[{i, j}, {i, 0, m - 1}, {j, 0, m - 1}], 1];
  dof = AssociationThread[nodes -> Range[0, m^2 - 1]];          (* node (0,0) -> 0: grounded *)
  els = Flatten[Table[Select[{{i + 1, j}, {i, j + 1}}, Max[#] < m &] /.
       nb : {_Integer, _Integer} :> {{dof[{i, j}], dof[nb]}, {{1, -1}, {-1, 1}}}, {i, 0, m - 1}, {j, 0, m - 1}], 2];
  F = ReplacePart[ConstantArray[0, m^2 - 1], dof[{1, 2}] -> 1];      (* unit current in at (1,2) *)
  {m^2 - 1, els, F}];

lShapeQ[a_, b_] := !(a >= 3 && b >= 3);
lShapeNodeQ[i_, j_] := And @@ Flatten[Table[0 <= a < 6 && 0 <= b < 6 && lShapeQ[a, b], {a, i - 1, i}, {b, j - 1, j}]];

problem[1] := bar[8, False];
problem[2] := bar[12, True];
problem[3] := q1[4, 4, 1/4, True &, interiorQ[4, 4]];
problem[4] := q1[8, 4, 1/4, True &, interiorQ[8, 4]];
problem[5] := q1[6, 6, 1/6, True &, interiorQ[6, 6], {3, 3}];
problem[6] := q1[6, 6, 1/6, True &, !(#2 == 0 && #1 <= 2) &];
problem[7] := gridGraph[4];
problem[8] := p1cell[4];
problem[9] := q1[6, 6, 1/6, lShapeQ, lShapeNodeQ];
problem[10] := q1[10, 10, 1/10, True &, interiorQ[10, 10]];
bounds = {1 -> {1, 10, 32}, 2 -> {10, 100, 606}, 3 -> {1, 10, 48}, 4 -> {1, 100, 1616},
   5 -> {1, 100, 1818}, 6 -> {1, 100, 360}, 7 -> {1, 10, 60}, 8 -> {1, 100, 808},
   9 -> {1, 100, 1350}, 10 -> {1, 100, 5050}};

(* the compliance F . U of the design cv (a numeric black box) *)
compliance[{nDof_, els_, F_}, cv_ /; VectorQ[cv, NumericQ]] := Module[{K, idx},
  K = ConstantArray[0 First[cv], {nDof, nDof}];
  Do[idx = els[[e, 1]];
     Do[If[idx[[r]] > 0 && idx[[s]] > 0,
          K[[idx[[r]], idx[[s]]]] += cv[[e]] els[[e, 2, r, s]]], {r, Length[idx]}, {s, Length[idx]}],
     {e, Length[els]}];
  F . LinearSolve[K, F]];

designMin[k_, opts___] := Module[{prob = problem[k], lo, hi, vol, m, cv},
  {lo, hi, vol} = k /. bounds;
  m = Length[prob[[2]]];
  cv = Array[cc, m];
  First[FindMinimum[{compliance[prob, cv], And @@ Thread[lo <= cv <= hi], Total[cv] == vol},
    Transpose[{cv, ConstantArray[vol/m, m]}], opts]]];

Do[With[{k = k},
    vrun["bd19." <> ToString[k] <> " machine", designMin[k], 10, 1, 90];
    vrun["bd19." <> ToString[k] <> " wp50", designMin[k, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 90]],
  {k, 1, 10}];
