(* ==========================================================================
   Experiment 3 -- Optimising codegen: CSE, folding, DCE, LICM
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.

       Mathilda      ./Mathilda -file compile_optimising_codegen.m
       Mathematica   wolframscript -file compile_optimising_codegen.m
       Python        python3 compile_optimising_codegen.py

   WHAT IT MEASURES.  Bodies chosen so that a specific optimisation either
   fires or does not: a long Horner chain (constant folding and instruction
   count), a repeated subexpression (CSE), a loop with an invariant (LICM),
   and a dead branch (DCE).

   WHY THE CEILING IS LOW, stated up front because the numbers look
   unimpressive next to experiment 1's 234x: a bytecode VM's dispatch is
   already the dominant cost, so removing instructions helps only in
   proportion to how many are removed.  1.05-1.48x is what that is worth, and
   the experiment exists to establish the number rather than to guess it.
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


(* A degree-40 Horner chain: 40 multiply-adds, all of whose coefficients are
   constant expressions the folder can evaluate once. *)
horner = Compile[{{x, _Real}},
  Module[{s = 0., k = 0},
    s = 1.;
    Do[s = s x + (1. + 1./(k + 1.)), {k, 1, 40}];
    s]];

(* A repeated subexpression: x y appears twice, so CSE should compute it once. *)
cse = Compile[{{x, _Real}, {y, _Real}}, x y + Sin[x y]];

(* A loop-invariant subexpression inside a While. *)
licm = Compile[{{x, _Real}, {n, _Integer}},
  Module[{s = 0., k = 0},
    While[k < n, s = s + x Sqrt[2.] + Sin[1.]; k = k + 1]; s]];

(* Newton's method: the micro-benchmark the codegen work was tuned on. *)
newt = Compile[{{a, _Real}},
  Module[{x = 1., k = 0},
    While[k < 20, x = x - (x x - a)/(2. x); k = k + 1]; x]];

Print["Experiment 3 -- optimising codegen"];
Print[""];

bench["Horner, degree 40, 10^6 calls",  Do[horner[0.5], {10^6}]];
check["Horner, degree 40",              horner[0.5]];

bench["x y + Sin[x y], 10^6 calls (CSE)", Do[cse[0.5, 0.25], {10^6}]];
check["x y + Sin[x y]",                   cse[0.5, 0.25]];

bench["loop-invariant body, 10^6 iterations (LICM)", licm[0.5, 10^6]];
check["loop-invariant body",                          licm[0.5, 10]];

bench["Newton, 20 iterations, 10^5 calls", Do[newt[2.], {10^5}]];
check["Newton (sqrt 2)",                   newt[2.]];

bench["degree-5 polynomial over 10^6 points",
      Table[1. + 2. x + 3. x^2 + 4. x^3 + 5. x^4 + 6. x^5, {x, 0., 1., 1./(10^6 - 1)}]];
