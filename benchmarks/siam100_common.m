(* ==========================================================================
   benchmarks/siam100_common.m -- shared by the SIAM 100-Digit Challenge
   experiments (166-195) and their holdouts (196-225).
   ==========================================================================

   Loaded after ../harness.m by every siam100 experiment:

       Get["../harness.m"]; Get["../siam100_common.m"];

   It reuses the ISSAC '97 helpers (sig, run, ruleVal, $cap) unchanged and adds
   one thing: the VALUE line.

   Every problem is ONE real number, and the question is how many of its digits
   each system gets.  The CHECK line (sig[v, d], d = the digits every correct
   system is expected to reach) is what run_all.py compares across systems, as
   in ISSAC.  The VALUE line carries every digit the result actually holds, and
   siam100_table.py scores it against the 60-digit reference in the
   experiment's truth.tsv:

       VALUE  <label>  {sign * leading digits, decimal exponent}

   run_all.py ignores VALUE (its LINE_RE matches only BENCH/CHECK/REQUIRE/SKIP),
   so the extra line costs the weekly job nothing.

   The case wrapper:

       vrun[label, expr, d]            one timed run
       vrun[label, expr, d, reps]      warm-up + min of reps
       vrun[label, expr, d, reps, cap] with a per-case cap in seconds

   expr must evaluate to the answer itself (a real number).  A FindRoot /
   FindMinimum / NMinimize result is unpacked by the case (ruleVal or First)
   inside expr, so that the timing includes nothing but the solve.
   ========================================================================== *)

Get["../issac97_common.m"];

(* The digits a value actually carries: 17 for a machine real (16 are honest,
   the 17th lets the scorer see the rounding), its Precision for an MPFR value,
   60 for an exact one (the scorer caps at the tier anyway). *)
availDigits[v_] := Which[
  MachineNumberQ[v], 17,
  Precision[v] === Infinity, 60,
  True, Max[1, Floor[Precision[v]]]];

(* A real answer, or a string meaning "not an answer".  A complex result with a
   nonzero imaginary part is not an answer to a real question; a complex
   result whose imaginary part is an exact or machine zero is its real part. *)
realOf[v_] := Which[
  StringQ[v], v,
  !NumericQ[v], "unevaluated",
  Head[v] === Complex && TrueQ[Im[v] == 0], Re[v],
  Head[v] === Complex, "complex",
  True, v];

SetAttributes[vrun, HoldRest];

vrun[label_String, expr_, d_Integer] := vrun[label, expr, d, 1, $cap];
vrun[label_String, expr_, d_Integer, reps_Integer] :=
  vrun[label, expr, d, reps, $cap];
vrun[label_String, expr_, d_Integer, reps_Integer, cap_] :=
  Module[{r, ts = {}, v},
    (* CheckAbort: an abort that escapes TimeConstrained (seen in
       wolframscript when NMaximize overruns its cap) must end this case
       only, not silence every later case in the file. *)
    If[reps > 1, CheckAbort[TimeConstrained[expr, cap, "TIMEOUT"], "TIMEOUT"]];
    Do[r = AbsoluteTiming[CheckAbort[TimeConstrained[expr, cap, "TIMEOUT"],
                                     "TIMEOUT"]];
       AppendTo[ts, First[r]], {reps}];
    v = If[Last[r] === "TIMEOUT", "timeout", realOf[Last[r]]];
    If[StringQ[v],
      check[label, v],
      Print["BENCH\t", label, "\t", ToString[Round[1000. Min[ts], 0.001]]];
      check[label, sig[v, d]];
      Print["VALUE\t", label, "\t",
            ToString[sig[v, availDigits[v]], InputForm]]]];
