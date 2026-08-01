(* ==========================================================================
   Experiment 17 -- Neural-network training: a chain of transposed GEMMs
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.

       Mathilda      ./Mathilda -file neural_network_training.m
       Mathematica   wolframscript -file neural_network_training.m
       Python        python3 neural_network_training.py

   WHAT IT MEASURES.  The third sweep's logistic-regression row has ONE weight
   matrix.  Real training has a chain, and the backward pass is a chain of
   TRANSPOSED products:

       dW2 = a1^T . d2      dh = d2 . W2^T      dW1 = X^T . d1

   Every one of those is a Transpose feeding a Dot, which is the operation the
   fourth sweep named as its largest remaining application gap.  A two-layer
   MLP is the smallest kernel that exercises the whole chain.

   BIASES ARE A COLUMN OF X, NOT A SEPARATE VECTOR.  That is what a
   performance-minded implementation does in all three languages, and it keeps
   the three columns identical: Wolfram's Plus threads the OUTER level, so
   `matrix + rowVector` is not a NumPy-style broadcast and the honest spelling
   would otherwise have to differ per system.

   DETERMINISM.  The timed run uses random data, which no two systems can
   share.  The check trains a DETERMINISTIC 256 x 17 -> 8 -> 4 network for 20
   steps and compares the cross-entropy loss, which depends on every weight
   and every gradient along the way.
   ========================================================================== *)

(* ---- shared reporting helpers (identical in every experiment file) ------- *)

SetAttributes[bench, HoldRest];

bench[label_String, expr_] := Module[{ts},
  expr;
  ts = Table[First[AbsoluteTiming[expr]], {3}];
  Print[StringPadRight[label, 52], ToString[Round[1000. Min[ts], 0.001]], " ms"]
];

check[label_String, value_] :=
  Print[StringPadRight[label, 52], "check = ", value];

relu[zz_] := Ramp[zz];

(* ---- the timed network -------------------------------------------------- *)

nnb  = 1024;      (* minibatch                                              *)
nni  = 785;       (* inputs, including the constant bias column             *)
nnh  = 128;       (* hidden units                                           *)
nno  = 10;        (* classes                                                *)
nnlr = 0.05;      (* learning rate                                          *)

nnX = Join[RandomReal[{-1, 1}, {nnb, nni - 1}], Table[{1.}, {nnb}], 2];
nnY = Table[If[j == Mod[i, nno] + 1, 1., 0.], {i, 1, nnb}, {j, 1, nno}];

(* Deterministic initial weights: no RNG, so the two CAS start identically. *)
nnW1 = Table[0.05 Sin[1. (i + nni j)], {i, 1, nni}, {j, 1, nnh}];
nnW2 = Table[0.05 Cos[1. (i + nnh j)], {i, 1, nnh}, {j, 1, nno}];

(* Row-wise softmax, shifted by the row maximum for numerical stability.

   Both threadings here rely on Wolfram's outer-level rule: `z - Map[Max, z]`
   pairs row i with the scalar max of row i, and `ex/Total[ex, {2}]` pairs row
   i with its own sum.  Neither is a NumPy broadcast; both are exactly what
   the outer-level rule gives. *)
sfmax[z_] := Module[{ex}, ex = Exp[z - Map[Max, z]]; ex/Total[ex, {2}]];

(* One SGD step: forward, softmax cross-entropy, backward, update. *)
mlpstep[w_] := Module[{z1, a1, z2, sm, d2, d1, g1, g2},
  z1 = nnX . w[[1]];                       (* 1024x785 . 785x128            *)
  a1 = relu[z1];
  z2 = a1 . w[[2]];                        (* 1024x128 . 128x10             *)
  sm = sfmax[z2];
  d2 = (sm - nnY)/nnb;                     (* dL/dz2 for cross-entropy      *)
  d1 = (d2 . Transpose[w[[2]]]) UnitStep[z1];   (* backprop through ReLU    *)
  g2 = Transpose[a1] . d2;
  g1 = Transpose[nnX] . d1;                (* Transpose[nnX] is LOOP-       *)
                                           (* INVARIANT and recomputed here *)
  {w[[1]] - nnlr g1, w[[2]] - nnlr g2}];

mlptrain[m_] := Nest[mlpstep, {nnW1, nnW2}, m];

(* ---- the deterministic check network ------------------------------------ *)

cb = 256; ci = 17; ch = 8; co = 4; clr = 0.1;
cX = Join[Table[Sin[0.3 (i + ci j)], {i, 1, cb}, {j, 1, ci - 1}],
          Table[{1.}, {cb}], 2];
cY = Table[If[j == Mod[i, co] + 1, 1., 0.], {i, 1, cb}, {j, 1, co}];
cW1 = Table[0.2 Sin[1. (i + ci j)], {i, 1, ci}, {j, 1, ch}];
cW2 = Table[0.2 Cos[1. (i + ch j)], {i, 1, ch}, {j, 1, co}];

cstep[w_] := Module[{z1, a1, z2, sm, d2, d1, g1, g2},
  z1 = cX . w[[1]]; a1 = relu[z1]; z2 = a1 . w[[2]];
  sm = sfmax[z2]; d2 = (sm - cY)/cb;
  d1 = (d2 . Transpose[w[[2]]]) UnitStep[z1];
  g2 = Transpose[a1] . d2; g1 = Transpose[cX] . d1;
  {w[[1]] - clr g1, w[[2]] - clr g2}];

closs[m_] := Module[{w, sm},
  w  = Nest[cstep, {cW1, cW2}, m];
  sm = sfmax[relu[cX . w[[1]]] . w[[2]]];
  -Total[cY Log[sm + 1.0*^-12], 2]/cb];

(* ---- inference ---------------------------------------------------------- *)

nnq  = 8192;
nnXq = Join[RandomReal[{-1, 1}, {nnq, nni - 1}], Table[{1.}, {nnq}], 2];
mlpinfer[] := Total[sfmax[relu[nnXq . nnW1] . nnW2], 2];

(* ---- run ---------------------------------------------------------------- *)

Print["Experiment 17 -- neural-network training"];
Print[""];

bench["MLP training, 785-128-10, batch 1024, 100 steps", mlptrain[100]];
check["MLP training (loss after 20 deterministic steps)", closs[20]];

bench["MLP inference, 8192 x 785 forward pass",          mlpinfer[]];
check["MLP inference (loss after 1 step)", closs[1]];

(* ---- one training step, split ------------------------------------------- *)
(*
   Two lines dominate, and neither is arithmetic: Transpose[nnX] is
   loop-invariant and recomputed every step, and the ReLU-derivative mask is a
   mixed float64 x int64 product.  Splitting the step is what shows that;
   see README.md for the roadmap it produces.
*)
Print[""];
Print["-- one training step, split (per step) --"];

sz1 = nnX . nnW1;
sa1 = relu[sz1];
sd2 = RandomReal[{0, 1}, {nnb, nno}];
sda = sd2 . Transpose[nnW2];
sdz = sda UnitStep[sz1];

bench["nnX . W1        (forward GEMM)",  nnX . nnW1];
bench["Ramp[z1]        (ReLU)",          relu[sz1]];
bench["Map[Max, z2]    (softmax shift)", Map[Max, sa1 . nnW2]];
bench["Transpose[nnX]  (LOOP-INVARIANT)", Transpose[nnX]];
bench["Transpose[a1]",                    Transpose[sa1]];
bench["Transpose[a1] . d2",               Transpose[sa1] . sd2];
bench["da1 UnitStep[z1] (mixed f64 x i64)", sda UnitStep[sz1]];
bench["Transpose[nnX] . dz1",             Transpose[nnX] . sdz];
bench["W1 - lr g1      (update)",         nnW1 - nnlr Transpose[nnX] . sdz];
