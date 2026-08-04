(* ==========================================================================
   Experiment 01 -- Symbolic integration of rational functions
   ==========================================================================

       cd benchmarks/01-integrate-rational
       ../../Mathilda -file integrate_rational.m
       wolframscript  -file integrate_rational.m
       python3 integrate_rational.py

   WHAT IT MEASURES.  The rational-function branch of `Integrate`: partial
   fraction decomposition, Hermite reduction for repeated denominator factors,
   and the Rothstein-Trager / Lazard-Rioboo-Trager step that produces the
   logarithmic part.  This is `src/calculus/intrat.c` and `src/parfrac.c`.

   WHY IT IS HERE.  Every one of the twenty experiments in docs/experiments/ is
   numeric, measured against NumPy.  The symbolic subsystems -- among the
   largest in the tree -- have never had an external baseline, so nothing in the
   repo could answer "is our Integrate slow?".  sympy is that baseline.

   HOW THE VALUE CHECK SURVIVES A DIFFERENT ANSWER.  Two systems can return
   algebraically different antiderivatives that are both correct: they differ by
   a constant, and by how the logarithmic part is grouped.  So the check is not
   the expression -- it is the DEFINITE integral recovered from it,
   F(b) - F(a), which is invariant under +C and under regrouping.  The interval
   is chosen to sit clear of every pole (ratFn[n] has poles at 1..n) so no
   branch cut is crossed.

   READ THIS ALONGSIDE ABSENT.md.  A symbolic row that is slow is usually slow
   because the other system has a *method* we do not, which is different work
   from a numeric fast path that was not taken.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"Integrate", "Apart", "Together", "D", "Factor"}];

(* ---- 1. Squarefree denominator: partial fractions + log part ------------ *)

f6 = ratFn[6, x];
bench["integrate rational, squarefree deg 6", Integrate[f6, x];];
check["integrate rational, squarefree deg 6",
  Round[10^6 ((Integrate[f6, x] /. x -> 7.5) - (Integrate[f6, x] /. x -> 6.5))]];

(* ---- 2. Repeated denominator factors: the Hermite reduction branch ------ *)

f4r = ratFnRepeated[4, x];
bench["integrate rational, repeated deg 4", Integrate[f4r, x];];
check["integrate rational, repeated deg 4",
  Round[10^6 ((Integrate[f4r, x] /. x -> 5.5) - (Integrate[f4r, x] /. x -> 4.5))]];

(* ---- 3. Irreducible denominator: no rational part, all logarithmic ------ *)

(* x^4 + 1 factors over Q into two irreducible quadratics, so the log part needs
   a genuine algebraic extension rather than a rational root, and the
   antiderivative is a closed form in Log and ArcTan.

   NOT x^5 + x + 1, which was the first choice: its antiderivative contains an
   unresolved root sum, and the two systems render those differently enough that
   substituting a number leaves `#1` slots behind on one side.  The TIMING would
   have been fine; the value check is what became unreliable, and a row that
   cannot be checked is a row that cannot be trusted. *)
fq = 1/(x^4 + 1);
bench["integrate rational, irreducible quartic", Integrate[fq, x];];
check["integrate rational, irreducible quartic",
  Round[10^6 N[(Integrate[fq, x] /. x -> 5/2) - (Integrate[fq, x] /. x -> 3/2)]]];

(* ---- 4. Polynomial integration: the trivial baseline -------------------- *)

(* Included deliberately as a floor.  If this row is slow, the cost is in
   expression handling rather than in any integration algorithm, which changes
   where the week's work goes. *)
p40 = denseUPoly[40, x];
bench["integrate polynomial, deg 40", Integrate[p40, x];];
check["integrate polynomial, deg 40",
  Round[10^6 N[Integrate[p40, x] /. x -> 1]]];

(* ---- 5. Partial fractions alone, without the integration ---------------- *)

(* Separates the decomposition cost from the integration cost: if row 1 is slow
   and this row is fast, the log part is the problem, not the decomposition. *)
f8 = ratFn[8, x];
bench["Apart, squarefree deg 8", Apart[f8, x];];
check["Apart, squarefree deg 8",
  Round[10^6 (Apart[f8, x] /. x -> 9.5)]];

(* ---- 6. Together: the inverse direction -------------------------------- *)

bench["Together, 8 summands", Together[Sum[1/(x - i), {i, 1, 8}]];];
check["Together, 8 summands",
  Round[10^6 (Together[Sum[1/(x - i), {i, 1, 8}]] /. x -> 9.5)]];
