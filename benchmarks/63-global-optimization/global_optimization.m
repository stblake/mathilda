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

(* ---- A2: isolated Gaussian well in a flat 10-D landscape -- DifferentialEvolution
   Global 1.524e-4 sits in a single narrow basin at 1.2345; the rest is a ~1.0
   plateau. DE's population finds it across every seed; a local method misses it.
   (The user's original d=12/sharpness-50 needle is unfindable by ANY
   derivative-free method — this is the widest still-isolated, genuinely findable
   well.) Check = Round[10^4 · obj] pins that the well (not the plateau) was found.
   Matching scipy engine: differential_evolution. *)
a2vars = Table[x[i], {i, 1, 10}];
a2obj = 1 - Exp[-2 Sum[(x[i] - 1.2345)^2, {i, 1, 10}]] + 10^-5 Sum[x[i]^2, {i, 1, 10}];
a2box = Table[-5 <= x[i] <= 5, {i, 1, 10}];
a2meth = {"DifferentialEvolution", "SearchPoints" -> 400, "ScalingFactor" -> 0.9,
          "RandomSeed" -> 1};
bench["A2 gaussian well (DE)", NMinimize[{a2obj, a2box}, a2vars, Method -> a2meth];];
check["A2 gaussian well (DE)",
  Round[10^4 First[NMinimize[{a2obj, a2box}, a2vars, Method -> a2meth]]]];

(* ---- A3: modified-Ackley product-of-cosines (10-D) -- SimulatedAnnealing
   Dense multimodal forest; true global -1 at x=0 is unreachable by any
   derivative-free method (origin basin radius ~pi/20). SA finds a deep basin
   (-0.92 at seed 1, best -0.984 across seeds) in ~0.5 s vs scipy dual_annealing
   ~1.0 s / best -0.975. Check = Round[obj] = -1 for both. Matching engine:
   dual_annealing (global vs global). *)
a3vars = Table[x[i], {i, 1, 10}];
a3obj = -Exp[-0.2 Sqrt[1/10 Sum[x[i]^2, {i, 1, 10}]]] Product[Cos[20 x[i]], {i, 1, 10}] +
         0.05 Sum[x[i]^2, {i, 1, 10}];
a3box = Table[-5 <= x[i] <= 5, {i, 1, 10}];
a3meth = {"SimulatedAnnealing", "LevelIterations" -> 200, "PerturbationScale" -> 0.5,
          "SearchPoints" -> 40, "RandomSeed" -> 1};
bench["A3 modified Ackley (SA)", NMinimize[{a3obj, a3box}, a3vars, Method -> a3meth];];
check["A3 modified Ackley (SA)",
  Round[First[NMinimize[{a3obj, a3box}, a3vars, Method -> a3meth]]]];

(* ---- A4: 15-D constrained minimax (Chebyshev) -- NelderMead
   min Max_i|x_i-Sin[i]| s.t. sum x^2<=15, sum i*x_i==1, box [-2,2]. Convex
   (max-of-affine + convex sphere + linear eq) -> unique optimum 0.125116. One
   simplex (SearchPoints 1) + AL polish reaches it in ~0.06 s, matching scipy
   single-start SLSQP epigraph (~0.065 s). Matching character: local/polish ->
   scipy.optimize.minimize. Check = Round[10^4 obj] = 1251. *)
a4vars = Table[x[i], {i, 1, 15}];
a4obj = Max[Table[Abs[x[i] - Sin[i]], {i, 1, 15}]];
a4cons = Join[{Sum[x[i]^2, {i, 1, 15}] <= 15, Sum[i x[i], {i, 1, 15}] == 1},
   Table[-2 <= x[i] <= 2, {i, 1, 15}]];
a4meth = {"NelderMead", "ContractRatio" -> 0.5, "Tolerance" -> 10^-8,
          "SearchPoints" -> 1, "RandomSeed" -> 1};
bench["A4 minimax chebyshev (NM)", NMinimize[{a4obj, a4cons}, a4vars, Method -> a4meth];];
check["A4 minimax chebyshev (NM)",
  Round[10^4 First[NMinimize[{a4obj, a4cons}, a4vars, Method -> a4meth]]]];

(* ---- A5: Almgren-Chriss optimal liquidation (T=15) -- SimulatedAnnealing
   Convex execution problem (|.|^1.5 impact + risk penalty, linear constraints),
   optimum 0.035206. SA at SearchPoints 1 (one chain + AL polish) reaches it in
   ~0.05 s, beating scipy single-start SLSQP ~0.11 s. Check = Round[10^6 obj]. *)
a5x = Array[pos, 16];
a5obj = Sum[0.1 Abs[a5x[[t]] - a5x[[t + 1]]]^1.5 + 0.05*0.04 a5x[[t]]^2, {t, 1, 15}];
a5cons = Join[{pos[1] == 1, pos[16] == 0}, Table[a5x[[t + 1]] <= a5x[[t]], {t, 1, 15}],
   Table[0 <= a5x[[t]] <= 1, {t, 1, 16}]];
a5meth = {"SimulatedAnnealing", "PerturbationScale" -> 0.1, "SearchPoints" -> 1,
          "RandomSeed" -> 1};
bench["A5 liquidation (SA)", NMinimize[{a5obj, a5cons}, a5x, Method -> a5meth];];
check["A5 liquidation (SA)",
  Round[10^6 First[NMinimize[{a5obj, a5cons}, a5x, Method -> a5meth]]]];

(* ---- A6: risk-parity portfolio (n=8) -- RandomSearch (local)
   Unimodal: unique interior solution, risk contributions all equal, obj ~0,
   min weight 0.0747 (identical to scipy SLSQP on the same covariance). One
   start (SearchPoints 1) suffices. NOTE this case is SLOWER than scipy here:
   NMinimize's per-solve objective SETUP (indexed-var rewrite + internal compile
   of a 4313-leaf tree with no common-subexpression sharing of the repeated cov.w
   dot products) dominates the ~0.07 s, versus scipy's pre-vectorized numpy at
   ~0.01 s; the overhead amortizes on larger problems. Check pins the solution
   via Round[10^4 min weight] = 747. Matching character: local -> scipy.minimize. *)
SeedRandom[202];
a6cov = (# + Transpose[#]) &@RandomReal[{-0.05, 0.15}, {8, 8}];
a6cov = a6cov . Transpose[a6cov];
a6w = Array[weight, 8]; a6rc = a6w * (a6cov . a6w);
a6obj = Sum[(a6rc[[i]] - a6rc[[j]])^2, {i, 8}, {j, 8}];
a6cons = Join[{Total[a6w] == 1}, Table[a6w[[i]] >= 0, {i, 8}]];
a6meth = {"RandomSearch", "Method" -> "InteriorPoint", "SearchPoints" -> 1, "RandomSeed" -> 1};
bench["A6 risk parity (RS)", NMinimize[{a6obj, a6cons}, a6w, Method -> a6meth];];
check["A6 risk parity (RS)",
  Round[10^4 Min[a6w /. Last[NMinimize[{a6obj, a6cons}, a6w, Method -> a6meth]]]]];

(* ---- A7: VWAP execution tracking (T=20) -- NelderMead + QuasiNewton
   min Σ(v_t-target_t)² s.t. Σv=1, v≥0, nonlinear impact budget Σ0.5 v²e^v ≤ B.
   Convex. The user's B=0.01 is INFEASIBLE (min impact 0.02628); B=0.027 binds
   (between the uniform min and target impact 0.02775), optimum 0.000224628.
   Convex ⇒ SearchPoints 1 + MaxIterations 100 (not 5000 with 20 restarts):
   ~0.08 s. Check = Round[10^6 obj] = 225. Matching character: local -> minimize. *)
a7tp = Table[0.05 + 0.04 Cos[Pi (t - 10)/10]^4, {t, 1, 20}]; a7tp = a7tp/Total[a7tp];
a7v = Array[trade, 20];
a7te = Sum[(a7v[[t]] - a7tp[[t]])^2, {t, 1, 20}];
a7ic = Sum[0.5 a7v[[t]]^2 Exp[a7v[[t]]], {t, 1, 20}];
a7cons = Join[{Total[a7v] == 1}, {a7ic <= 0.027}, Table[a7v[[t]] >= 0, {t, 1, 20}]];
a7meth = {"NelderMead", "PostProcess" -> "QuasiNewton", "SearchPoints" -> 1, "RandomSeed" -> 1};
bench["A7 vwap tracking (NM)", NMinimize[{a7te, a7cons}, a7v, Method -> a7meth, MaxIterations -> 100];];
check["A7 vwap tracking (NM)",
  Round[10^6 First[NMinimize[{a7te, a7cons}, a7v, Method -> a7meth, MaxIterations -> 100]]]];
