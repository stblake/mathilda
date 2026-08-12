(* Experiment 24 -- The heterogeneous return.
   ROADMAP ITEM 1 -- "a heterogeneous packed tuple: a container holding n buffers
   of independent shape and dtype without materialising them", valued at 4.07x on
   MLP training and 120x on a mixed-dtype return (experiments 13 and 17).

   THE DEFECT.  A function returning {reals, ints, mask} destroys all three
   arrays: one is an integer buffer, one is boolean, and a mixed-dtype list of
   packed rows cannot be absorbed, so every element is materialised as a boxed
   Expr -- 96x at the return and 120x on the caller's NEXT operation.

   Each row below returns a differently-typed tuple and then CONSUMES it, because
   the second cost only appears on the consumer. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Total", "Length", "Part", "Map", "UnitStep", "Round", "Floor"}];

n = 500000;
v = rand01[{n}];

(* 1. Homogeneous return: the control.  All-real, should stay packed. *)
homog[x_] := {x + 1., x * 2.};
bench["return {real, real}, then Total",
  Module[{p = homog[v]}, Total[p[[1]]] + Total[p[[2]]]];];
check["return {real, real}, then Total",
  Round[10^4 N[Total[First[homog[Range[10.]]]]]]];

(* 2. Mixed real + integer: the defect. *)
mixed2[x_] := {x + 1., Round[x * 10]};
bench["return {real, int}, then Total",
  Module[{p = mixed2[v]}, Total[p[[1]]] + Total[p[[2]]]];];
check["return {real, int}, then Total",
  Round[10^4 N[Total[Last[mixed2[Range[10.]/10]]]]]];

(* 3. Real + integer + boolean mask: the exact shape from experiment 13. *)
mixed3[x_] := {x + 1., Round[x * 10], UnitStep[x - 1/2]};
bench["return {real, int, mask}, then Total",
  Module[{p = mixed3[v]},
    Total[p[[1]]] + Total[p[[2]]] + Total[p[[3]]]];];
check["return {real, int, mask}, then Total",
  Round[10^4 N[Total[Last[mixed3[Range[10.]/10]]]]]];

(* 4. The return alone, without consuming it: isolates the two costs. *)
bench["return {real, int, mask}, discarded", mixed3[v];];
check["return {real, int, mask}, discarded", Length[mixed3[Range[10.]/10]]];

(* 5. Ragged shapes as well as mixed dtypes. *)
ragged[x_] := {x, Take[x, 1000], Round[Take[x, 100] * 10]};
bench["return ragged {n, 1000, 100}, then Total",
  Module[{p = ragged[v]},
    Total[p[[1]]] + Total[p[[2]]] + Total[p[[3]]]];];
check["return ragged {n, 1000, 100}, then Total",
  Length[ragged[rand01[{2000}]]]];
