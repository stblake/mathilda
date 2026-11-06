(* ==========================================================================
   benchmarks/issac97_common.m -- shared by experiments 126-135, the ISSAC '97
   Systems Challenge (one experiment per problem).
   ==========================================================================

   Loaded after ../harness.m by every issac97 experiment:

       Get["../harness.m"]; Get["../issac97_common.m"];

   The problems' answers are huge or tiny (a condition number of 1.77e389, an
   integral of 1.10e36300, a coefficient of 3.97e1426), so a check is the
   leading SIGNIFICANT digits and the decimal exponent, never Round[10^6 v]:

       sig[v, d]  ->  {sign * leading d digits as an integer, exponent}

   e.g. sig[21.19324037771, 8] = {21193240, 1}.  The digit count d is chosen per
   case at or below the precision the method guarantees, so the two systems
   agree whenever both are right, and a wrong answer cannot hide.

   A case that cannot finish is capped with TimeConstrained so that one runaway
   case (a 270 s Limit, a > 11 min NSum at 40 digits) costs
   at most $cap seconds instead of the whole file.  A timed-out or unevaluated
   case emits ONLY its CHECK line ("timeout" / "unevaluated"); with no BENCH
   line the runner classifies it INCOMPLETE for that system and never pairs a
   timing with a non-answer.
   ========================================================================== *)

$cap = 120;

(* Leading significant digits + exponent; strings for the non-answers. *)
sig[v_, d_Integer] := Which[
  v === "TIMEOUT", "timeout",
  !NumericQ[v], "unevaluated",
  (* TrueQ: a number with no significant digits (an accuracy-only zero, or a
     result whose precision collapsed) leaves v == 0 undecided, and a bare
     Which would then return itself as the "check value". *)
  TrueQ[v == 0] || Precision[v] === 0. || Precision[v] === 0, {0, 0},
  (* Never N[v]: on an arbitrary-precision number N[v] drops to machine
     precision, and every check past the 16th digit would compare rounding
     noise.  Exact input gets d + 10 digits; inexact input keeps its own. *)
  True, Module[{a = If[Precision[v] === Infinity, N[Abs[v], d + 10], Abs[v]], e},
    e = Floor[Log10[a]];
    {Sign[v] Round[a/10^(e - d + 1)], e}]];

(* run[label, expr, chk]  -- one timed run, no warm-up (every ISSAC case costs
   at least milliseconds of real work, most of them seconds).
   run[label, expr, chk, reps] -- warm-up + min of reps, for the fast rows.
   chk is a function of the result; it returns either the value to CHECK or
   a string ("unevaluated", "timeout", "failed") meaning "not an answer". *)

SetAttributes[run, HoldRest];

run[label_String, expr_, chk_] := run[label, expr, chk, 1, $cap];
run[label_String, expr_, chk_, reps_Integer] := run[label, expr, chk, reps, $cap];
run[label_String, expr_, chk_, reps_Integer, cap_] :=
  Module[{r, ts = {}, c},
    If[reps > 1, TimeConstrained[expr, cap, "TIMEOUT"]];        (* warm-up *)
    Do[r = AbsoluteTiming[TimeConstrained[expr, cap, "TIMEOUT"]];
       AppendTo[ts, First[r]], {reps}];
    c = If[Last[r] === "TIMEOUT", "timeout", chk[Last[r]]];
    If[!StringQ[c],
       Print["BENCH\t", label, "\t", ToString[Round[1000. Min[ts], 0.001]]]];
    check[label, c]];

(* A FindRoot / FindMinimum result -> the value of the variable, or "failed". *)
ruleVal[r_, x_] := If[MatchQ[r, {___Rule} | {_?NumericQ, {___Rule}}],
  x /. If[MatchQ[r, {_, {___Rule}}], Last[r], r], "failed"];
