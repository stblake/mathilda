(* --------------------------------------------------------------------
 * distributions.m -- HeavisideTheta / DiracDelta evaluation rules.
 *
 * These are the distributional functions produced by DSolve's Green's-
 * function / impulse solutions of forced linear ODEs.  The sifting property
 * of DiracDelta under a definite integral lives in C (src/calculus/
 * integrate_dirac.c) and the derivative rule HeavisideTheta' = DiracDelta is
 * in src/calculus/deriv.c; here we only add the value rules that fire on a
 * definite-sign NUMERIC argument, so that fitting an initial condition at the
 * base point resolves a shifted HeavisideTheta[t - a] / DiracDelta[t - a].
 * A symbolic argument is left inert.
 *
 * HeavisideTheta is taken left-continuous, H(0) = 0, so that a causal impulse
 * response x(t) = G(t, a) HeavisideTheta[t - a] satisfies its pre-impulse
 * initial conditions x(0) = x'(0) = 0 (matching the DSolve output for e.g.
 * x'' + 4 x == DiracDelta[t]).
 * -------------------------------------------------------------------- *)

HeavisideTheta[0] = 0;
HeavisideTheta[x_?Positive] := 1;
HeavisideTheta[x_?Negative] := 0;

DiracDelta[x_?Positive] := 0;
DiracDelta[x_?Negative] := 0;

(* The multivariate forms are products: HeavisideTheta[x1, x2, ...] is 1 when
 * every xi is positive and 0 when any is negative (or 0, by the convention
 * above); DiracDelta[x1, x2, ...] is 0 as soon as one xi is nonzero. *)
HeavisideTheta[x_, y__] /; AnyTrue[{x, y}, TrueQ[Negative[#]] || TrueQ[# == 0] &] := 0;
HeavisideTheta[x_, y__] /; AllTrue[{x, y}, TrueQ[Positive[#]] &] := 1;
HeavisideTheta[a___, p_?Positive, b___] /; Length[{a, b}] > 0 := HeavisideTheta[a, b];
DiracDelta[x_, y__] /; AnyTrue[{x, y}, TrueQ[Positive[#]] || TrueQ[Negative[#]] &] := 0;

(* Attributes as in the reference: Listable (FeynCalc's FCLoopFindSectors
 * reads sectors as HeavisideTheta[pows - 1/2] over a list), Orderless,
 * Protected.  Set last, so the rules above are installed first. *)
SetAttributes[{HeavisideTheta, DiracDelta}, {Listable, Orderless, Protected}];
