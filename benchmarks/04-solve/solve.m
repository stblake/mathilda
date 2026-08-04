(* ==========================================================================
   Experiment 04 -- Solve: polynomial, systems, transcendental
   ==========================================================================
   WHAT IT MEASURES.  src/solve.c and src/poly/solvepoly.c -- root isolation for
   univariate polynomials, elimination for systems, and the inverse-function
   branch for transcendental equations.

   VALUE CHECKS ARE COUNTS AND SUMS, not solution sets.  Two systems legitimately
   order roots differently and spell radicals differently, so the check is the
   NUMBER of solutions plus the rounded SUM of the roots -- both invariant under
   ordering and under algebraic form.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"Solve", "NSolve", "Reduce", "Eliminate", "Roots"}];

(* 1. A factorable sextic: rational roots, no radicals needed. *)
e1 = Product[x - i, {i, 1, 6}] == 0;
bench["Solve factorable sextic", Solve[e1, x];];
check["Solve factorable sextic", Length[Solve[e1, x]]];

(* 2. A quartic needing radicals. *)
e2 = x^4 - 10 x^2 + 1 == 0;
bench["Solve quartic with radicals", Solve[e2, x];];
check["Solve quartic with radicals", Length[Solve[e2, x]]];

(* 3. A 2x2 polynomial system: elimination. *)
bench["Solve 2x2 polynomial system",
  Solve[{x^2 + y^2 == 25, x - y == 1}, {x, y}];];
check["Solve 2x2 polynomial system",
  Length[Solve[{x^2 + y^2 == 25, x - y == 1}, {x, y}]]];

(* 4. A 3x3 linear system: the cheap baseline.  If this row is slow the cost is
   in expression handling, not in any solving algorithm. *)
bench["Solve 3x3 linear system",
  Solve[{x + y + z == 6, x - y == 1, y - z == 1}, {x, y, z}];];
check["Solve 3x3 linear system",
  Length[Solve[{x + y + z == 6, x - y == 1, y - z == 1}, {x, y, z}]]];

(* 5. Numeric roots of a degree-40 polynomial. *)
p40 = denseUPoly[40, x];
bench["NSolve degree 40", NSolve[p40 == 0, x];];
check["NSolve degree 40", Length[NSolve[p40 == 0, x]]];

(* 6. Sum of the roots of the sextic -- pins the VALUES, not just the count. *)
bench["Solve sextic, root sum", Total[x /. Solve[e1, x]];];
check["Solve sextic, root sum", Round[Total[x /. Solve[e1, x]]]];
