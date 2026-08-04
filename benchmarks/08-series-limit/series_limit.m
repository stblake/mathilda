(* ==========================================================================
   Experiment 08 -- Series expansion and limits
   ==========================================================================
   WHAT IT MEASURES.  src/calculus/series.c and the Gruntz limit machinery
   (src/calculus/gruntz.c, limit.c, limit_osc.c).

   Series and Limit fail in opposite directions, which is why both are here.
   Series is a cost problem: the answer is never in doubt, only how long the
   coefficient recursion takes.  Limit is a correctness problem: the hard cases
   need a comparability ordering, and a system without one returns the input
   unevaluated -- fast and wrong.

   Checks are single COEFFICIENTS and limit VALUES, both exact.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"Series", "SeriesCoefficient", "Limit", "Normal", "O", "Direction"}];

(* 1. A composed transcendental to order 20: coefficient recursion. *)
bench["Series Exp[Sin[x]] to order 20", Series[Exp[Sin[x]], {x, 0, 20}];];
check["Series Exp[Sin[x]] to order 20",
  SeriesCoefficient[Series[Exp[Sin[x]], {x, 0, 20}], 12]];

(* 2. A rational generating function to order 60: Fibonacci coefficients. *)
bench["Series 1/(1-x-x^2) to order 60", Series[1/(1 - x - x^2), {x, 0, 60}];];
check["Series 1/(1-x-x^2) to order 60",
  SeriesCoefficient[Series[1/(1 - x - x^2), {x, 0, 60}], 50]];

(* 3. Log of a series: needs the log of a unit power series. *)
bench["Series Log[1+Sin[x]] to order 24",
  Series[Log[1 + Sin[x]], {x, 0, 24}];];
check["Series Log[1+Sin[x]] to order 24",
  SeriesCoefficient[Series[Log[1 + Sin[x]], {x, 0, 24}], 15]];

(* 4. The textbook limit -- the floor for this subsystem. *)
bench["Limit Sin[x]/x at 0", Limit[Sin[x]/x, x -> 0];];
check["Limit Sin[x]/x at 0", Limit[Sin[x]/x, x -> 0]];

(* 5. A limit needing repeated cancellation of competing rates. *)
bench["Limit (Exp[x]-1-x)/x^2 at 0",
  Limit[(Exp[x] - 1 - x)/x^2, x -> 0];];
check["Limit (Exp[x]-1-x)/x^2 at 0", Limit[(Exp[x] - 1 - x)/x^2, x -> 0]];

(* 6. A Gruntz-class limit at infinity: exp and log at comparable rates.  This
   is the row that separates a real comparability ordering from a lucky guess. *)
bench["Limit x (Log[x+1]-Log[x]) at Infinity",
  Limit[x (Log[x + 1] - Log[x]), x -> Infinity];];
check["Limit x (Log[x+1]-Log[x]) at Infinity",
  Limit[x (Log[x + 1] - Log[x]), x -> Infinity]];
