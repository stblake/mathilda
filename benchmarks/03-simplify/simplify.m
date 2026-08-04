(* ==========================================================================
   Experiment 03 -- Simplify and FullSimplify
   ==========================================================================
   WHAT IT MEASURES.  src/simp/ -- the search-based simplifier: trig canonical
   forms, radical denesting, log/exp collapse, and rational cancellation.

   WHY SIMPLIFY IS THE HARDEST THING HERE TO BENCHMARK HONESTLY.  Simplify has
   no unique right answer, only a smaller one, so two systems can both "succeed"
   and return different expressions of different sizes.  A timing alone would
   reward whichever system gives up soonest.  So every row here checks a
   NUMERIC value at a rational point: any correct simplification preserves it,
   and a system that returned the input unchanged still passes the check while
   showing its cost -- which is the honest comparison.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"Simplify", "FullSimplify", "TrigReduce", "TrigExpand", "Cancel",
         "PowerExpand", "Factor", "Expand"}];

at[e_, v_] := Round[10^6 N[e /. x -> v]];

(* 1. A trig identity that collapses to 1. *)
t1 = Sin[x]^2 + Cos[x]^2;
bench["Simplify trig Pythagorean", Simplify[t1];];
check["Simplify trig Pythagorean", at[Simplify[t1], 1/3]];

(* 2. A multiple-angle expansion that must contract. *)
t2 = 8 Sin[x]^4 - 8 Sin[x]^2 + 1;
bench["Simplify quartic-to-Cos[4x]", Simplify[t2];];
check["Simplify quartic-to-Cos[4x]", at[Simplify[t2], 1/3]];

(* 3. Rational cancellation: a removable singularity. *)
t3 = (x^3 - 1)/(x - 1);
bench["Cancel removable singularity", Cancel[t3];];
check["Cancel removable singularity", at[Cancel[t3], 5/2]];

(* 4. A nested radical that denests. *)
t4 = Sqrt[3 + 2 Sqrt[2]];
bench["FullSimplify nested radical", FullSimplify[t4];];
check["FullSimplify nested radical", Round[10^6 N[FullSimplify[t4]]]];

(* 5. Log/exp collapse. *)
t5 = Log[Exp[x]] + Exp[Log[x]];
bench["Simplify log-exp collapse", Simplify[t5];];
check["Simplify log-exp collapse", at[Simplify[t5], 7/4]];

(* 6. TrigReduce on a product of four sines: linearisation. *)
t6 = Sin[x] Sin[2 x] Sin[3 x] Sin[4 x];
bench["TrigReduce product of 4 sines", TrigReduce[t6];];
check["TrigReduce product of 4 sines", at[TrigReduce[t6], 1/3]];
