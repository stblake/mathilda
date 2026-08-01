(* ==========================================================================
   Experiment 2 -- Compiled machine arrays and fused elementwise loops
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.

       Mathilda      ./Mathilda -file compile_array_fusion.m
       Mathematica   wolframscript -file compile_array_fusion.m
       Python        python3 compile_array_fusion.py

   WHAT IT MEASURES.  An elementwise expression over a 10^6-element vector,
   written three ways: interpreted, compiled, and compiled with the loop
   FUSED.  Fusion is the whole subject -- `v^2 + 2 v + 1` unfused is four
   passes over 8 MB with three temporaries; fused it is one pass with none.

   The Sin/Exp/Sqrt row is deliberately included because it is libm-bound
   rather than memory-bound, so it shows what fusion is worth when the
   arithmetic, not the traffic, is the cost -- less, and the file says so.
   ========================================================================== *)

SetAttributes[bench, HoldRest];

(* bench[label, expr] -- one untimed warm-up, then the MINIMUM of three timed
   runs.  The minimum, not the mean: we are measuring the cost of the work, and
   every source of noise on a loaded machine can only add. *)
bench[label_String, expr_] := Module[{ts},
  expr;
  ts = Table[First[AbsoluteTiming[expr]], {3}];
  Print[StringPadRight[label, 52], ToString[Round[1000. Min[ts], 0.001]], " ms"]
];

check[label_String, value_] :=
  Print[StringPadRight[label, 52], "check = ", value];


n = 10^6;
v = RandomReal[{0, 1}, n];

(* Compiled, fused: one pass, no temporaries. *)
f1 = Compile[{{a, _Real, 1}}, a^2 + 2 a + 1];
f2 = Compile[{{a, _Real, 1}}, Total[a^2 + 2 a + 1]];
f3 = Compile[{{a, _Real, 1}}, Total[a + a a]];
f4 = Compile[{{a, _Real, 1}}, Total[Sin[a] Exp[-a] + Sqrt[a]]];
f5 = Compile[{{a, _Real, 1}}, Sqrt[a] + a^2];
f6 = Compile[{{a, _Real, 1}}, Sin[a] Exp[-a] + Sqrt[a]];

Print["Experiment 2 -- compiled machine arrays and fused elementwise loops"];
Print[""];
Print["-- the same expression, interpreted then compiled --"];

bench["v^2 + 2 v + 1            (interpreted)", v^2 + 2 v + 1];
bench["v^2 + 2 v + 1            (compiled, fused)", f1[v]];
bench["Total[v^2 + 2 v + 1]     (interpreted)", Total[v^2 + 2 v + 1]];
bench["Total[v^2 + 2 v + 1]     (compiled, fused)", f2[v]];
bench["Total[v + v v]           (compiled, fused)", f3[v]];
bench["Sqrt[v] + v^2            (compiled, fused)", f5[v]];

Print[""];
Print["-- libm-bound bodies: fusion is worth less, because the cost is not traffic --"];
bench["Sin[v] Exp[-v] + Sqrt[v] (interpreted)", Sin[v] Exp[-v] + Sqrt[v]];
bench["Sin[v] Exp[-v] + Sqrt[v] (compiled, fused)", f6[v]];
bench["Total[Sin[v] Exp[-v] + Sqrt[v]] (compiled)", f4[v]];

Print[""];
Print["-- the elementwise primitives, for scale --"];
bench["Sin[v]",  Sin[v]];
bench["Exp[v]",  Exp[v]];
bench["v + 3. v", v + 3. v];

check["Total[v^2 + 2 v + 1] agrees", Round[f2[v] - Total[v^2 + 2 v + 1], 1.*^-6]];
