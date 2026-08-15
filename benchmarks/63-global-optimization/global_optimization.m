(* Experiment 63 -- Global / constrained optimization (NMinimize).

   MEASURES src/numerical_calculus/findmin.c: the derivative-free global engines
   (DifferentialEvolution / SimulatedAnnealing / NelderMead / RandomSearch) and
   the augmented-Lagrangian constrained local polish, against scipy's global
   optimizers (differential_evolution / dual_annealing) and, where the method
   character is local, scipy.optimize.minimize. Each case is timed in both
   systems and joined on its label; the check is the objective value at the
   optimum, rounded, which is invariant under which feasible optimum was found.

   Random-data problems bake IDENTICAL numeric data into both halves (Mathilda's
   RNG != numpy's, so the seed cannot be shared -- only the data). Cases are
   added one at a time as the optimization testbed works through them. *)

Get["../harness.m"];
Get["../data.m"];

require[{"NMinimize", "NMaximize"}];

(* ---- A1: refinery pooling (bilinear equality) -- DifferentialEvolution ----
   8 vars, two equality constraints (one BILINEAR: 3 cA + cB = px xA + py yA),
   linear capacity inequalities, box bounds. Global optimum -3900. The bilinear
   equality is what the augmented-Lagrangian polish upgrade was built for; DE
   reaches the global optimum across every seed. Matching scipy engine:
   differential_evolution (global vs global). *)
a1cost = 6 cA + 16 cB + 10 (xA + xB) - 9 (xA + yA) - 15 (xB + yB);
a1cons = {cA + cB - (px + py) == 0, 3 cA + 1 cB - px xA - py yA == 0,
   xA + yA <= 100, xB + yB <= 200, cA <= 300, cB <= 300, px <= 2.5, py <= 1.5,
   cA >= 0, cB >= 0, xA >= 0, xB >= 0, yA >= 0, yB >= 0, px >= 0, py >= 0,
   xA <= 100, yA <= 100, xB <= 200, yB <= 200};
a1vars = {px, py, xA, xB, yA, yB, cA, cB};
a1meth = {"DifferentialEvolution", "RandomSeed" -> 1};
bench["A1 refinery pooling (DE)", NMinimize[{a1cost, a1cons}, a1vars, Method -> a1meth];];
check["A1 refinery pooling (DE)",
  Round[First[NMinimize[{a1cost, a1cons}, a1vars, Method -> a1meth]]]];
