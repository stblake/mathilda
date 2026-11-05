(* ==========================================================================
   Experiment 123 -- native kernels against the Mathilda code they replace
   ==========================================================================

       cd benchmarks/123-native-vs-m
       ../../Mathilda -file native_vs_m.m
       wolframscript  -file native_vs_m.m

   WHAT IT MEASURES.  Each pair here is ONE computation written twice: as
   ordinary Wolfram-language code, which runs through the evaluator, and as a
   native C builtin (an Experimental` head) that must return the identical
   value.  The ".m" rows run in every system that reads this file.  The "C" rows
   run only where the head exists, and elsewhere print SKIP.  The ratio between
   the two rows on the same system is what the C kernel buys over the
   interpreter.  The ratio of the ".m" rows across systems compares evaluators.

   A new pair is added by appending a section in the same shape:
   a reference definition, bench rows for both, and a check per size.

   --------------------------------------------------------------------------
   PAIR 1: the recursive Dirac trace (Rolf Mertig, Woxi issue #2142)
   --------------------------------------------------------------------------
   diracTr below is the issue's definition verbatim.  Its native twin is
   Experimental`FeynCalc`DiracTr (src/experimental/feyncalc/diractr.c), which is
   SameQ to it, term order included.  n = 14 gives 13!! = 135135 terms in 91
   distinct g.

   Reported in the issue, all runs on one Linux machine, n = 14:
       Woxi v0.3.0        trace 40.82 s   Cases 1.27 s   1704 MB
       wolframscript 15   trace  3.42 s   Cases 0.30 s    100 MB
       Mathilda 0.6177    trace  1.52 s   Cases 0.19 s    101 MB
       FORM 5.0.2         0.58 s (same algorithm, medium.frm, 8.6 MB used)
   run_all.py has no FORM runner.  medium.frm sits beside this file so the
   number can be reproduced by hand with  time form -q medium.frm.
   ========================================================================== *)

Get["../harness.m"];

require[{"Experimental`FeynCalc`DiracTr"}];

(* Decide presence HERE, before any statement below is parsed.  Parsing a
   statement that spells Experimental`FeynCalc`DiracTr interns the symbol, so
   harness benchIf's run-time Names[] test would say "present" in a system
   without the builtin, and time an unevaluated call.  (Nor can benchIf take
   the name from a variable: it is HoldRest, and the held symbol matches no
   head_String, so the call is silently left unevaluated.) *)
$haveNative = haveQ["Experimental`FeynCalc`DiracTr"];
SetAttributes[{benchC, checkC}, HoldRest];
benchC[label_String, expr_] :=
  If[$haveNative, bench[label, expr], Print["SKIP\t", label, "\tExperimental`FeynCalc`DiracTr"]];
checkC[label_String, value_] := If[$haveNative, check[label, value], Null];

diracTr[] := 4;
diracTr[a_, b_] := 4 g[a, b];
diracTr[a_, rest___] := Module[{list = {rest}, n},
  Sum[(-1)^(m + 1) * Expand[g[a, list[[m]]] * diracTr @@ Delete[list, m]],
    {m, 1, Length[list]}]];

(* One scalar that pins the whole sum: its length, the distinct g it uses, and
   the sum of its coefficients. *)
dtDigest[r_] := {Length[r], Length[Union[Cases[r, _g, Infinity]]], r /. g[_, _] -> 1};

Do[
  With[{k = k},
    bench["diracTr .m n=" <> ToString[k], rM = diracTr @@ Array[x, k]];
    check["diracTr .m n=" <> ToString[k], dtDigest[rM]];
    benchC["DiracTr C n=" <> ToString[k],
      rC = Experimental`FeynCalc`DiracTr @@ Array[x, k]];
    checkC["DiracTr C n=" <> ToString[k], dtDigest[rC]];
    checkC["DiracTr C === .m n=" <> ToString[k], rC === rM];
  ],
  {k, {8, 10, 12, 14}}];

(* The issue's second timing, on the n = 14 results of both. *)
bench["Cases _g on .m n=14", Length[Union[Cases[rM, _g, Infinity]]]];
check["Cases _g on .m n=14", Length[Union[Cases[rM, _g, Infinity]]]];
benchC["Cases _g on C n=14", Length[Union[Cases[rC, _g, Infinity]]]];

(* Past what the recursion runs in benchmark time: the native kernel alone. *)
benchC["DiracTr C n=16", rC = Experimental`FeynCalc`DiracTr @@ Array[x, 16]];
checkC["DiracTr C n=16", dtDigest[rC]];

(* Repeated indices: equal factors fold to powers and equal terms combine. *)
rep = {mu, nu, rho, sig, mu, nu, rho, sig, a, b, a, b};
bench["diracTr .m repeated n=12", rM = diracTr @@ rep];
check["diracTr .m repeated n=12", {Length[rM], rM /. g[_, _] -> 1}];
benchC["DiracTr C repeated n=12", rC = Experimental`FeynCalc`DiracTr @@ rep];
checkC["DiracTr C repeated n=12", {Length[rC], rC /. g[_, _] -> 1}];
checkC["DiracTr C === .m repeated n=12", rC === rM];
