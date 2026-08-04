(* ==========================================================================
   Experiment 02 -- Symbolic integration, transcendental cases
   ==========================================================================

       cd benchmarks/02-integrate-transcendental
       ../../Mathilda -file integrate_transcendental.m

   WHAT IT MEASURES.  The Risch machinery: building a differential field tower
   over exp / log / tan, deciding elementary integrability, and returning a
   special function when the answer is not elementary.  This is the large
   `src/calculus/risch_*.c` family plus `integrate_risch_transcendental.c`.

   WHY THESE SIX.  Each one exercises a different branch, so a slow row names a
   subsystem rather than "integration":

     exp tower       x Exp[x^2]           substitution, no tower needed
     log tower       Log[x]^3             repeated integration by parts
     mixed tower     Sin[x] Exp[x]        two extensions, coupled
     tan tower       Tan[x]^3             the hyperbolic-tangent frontend
     rational-in-exp 1/(1 + Exp[x])       a rational function OF an exponential
     non-elementary  Exp[x]/x             must return ExpIntegralEi, not fail

   THE LAST ROW IS THE INTERESTING ONE.  Exp[x]/x has no elementary
   antiderivative.  A correct system proves that and returns Ei; a system
   without the proof either returns the integral unevaluated (fast, and wrong)
   or searches forever.  The value check is what tells those apart -- an
   unevaluated return costs almost nothing and would otherwise read as a win.

   VALUE CHECKS.  F(b) - F(a) recovered from the antiderivative, which is
   invariant under +C and under regrouping, evaluated away from any singularity.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"Integrate", "Exp", "Log", "Sin", "Tan", "ExpIntegralEi", "D"}];

dif[F_, a_, b_] := Round[10^6 N[(F /. x -> b) - (F /. x -> a)]];

(* ---- 1. exp tower: a plain substitution -------------------------------- *)

bench["integrate x Exp[x^2]", Integrate[x Exp[x^2], x];];
check["integrate x Exp[x^2]", dif[Integrate[x Exp[x^2], x], 1/4, 3/4]];

(* ---- 2. log tower: repeated integration by parts ----------------------- *)

bench["integrate Log[x]^3", Integrate[Log[x]^3, x];];
check["integrate Log[x]^3", dif[Integrate[Log[x]^3, x], 3/2, 5/2]];

(* ---- 3. mixed tower: two coupled extensions ---------------------------- *)

bench["integrate Sin[x] Exp[x]", Integrate[Sin[x] Exp[x], x];];
check["integrate Sin[x] Exp[x]", dif[Integrate[Sin[x] Exp[x], x], 1/2, 3/2]];

(* ---- 4. tan tower ------------------------------------------------------ *)

bench["integrate Tan[x]^3", Integrate[Tan[x]^3, x];];
check["integrate Tan[x]^3", dif[Integrate[Tan[x]^3, x], 1/4, 1/2]];

(* ---- 5. a rational function of an exponential -------------------------- *)

bench["integrate 1/(1+Exp[x])", Integrate[1/(1 + Exp[x]), x];];
check["integrate 1/(1+Exp[x])", dif[Integrate[1/(1 + Exp[x]), x], 1/2, 3/2]];

(* ---- 6. non-elementary: must return Ei, not give up -------------------- *)

(* If this returns unevaluated it is FAST and WRONG.  The check is the only
   thing standing between that and a reported win. *)
bench["integrate Exp[x]/x (non-elementary)", Integrate[Exp[x]/x, x];];
check["integrate Exp[x]/x (non-elementary)",
  dif[Integrate[Exp[x]/x, x], 3/2, 5/2]];
