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

(* DiracDelta is even, and the reference writes an argument whose terms are
 * mostly negated with the sign taken out: DiracDelta[1 - x1 - x2] is
 * DiracDelta[-1 + x1 + x2] and DiracDelta[-Sin[x]] is DiracDelta[Sin[x]],
 * while DiracDelta[a - b], DiracDelta[-a + b] and DiracDelta[1 - 2 x] stay.
 * The count: a term -t (coefficient exactly -1) against a term t (coefficient
 * exactly 1), numbers aside; more of the former flips.  An inexact argument is
 * left as written (DiracDelta[1.5 - x] stays).  FeynCalc's
 * FCFeynmanParameterJoin writes its delta function this way. *)
DiracDelta`Private`flipQ[e_] := FreeQ[e, _Real | _Complex] && Module[{t = If[Head[e] === Plus, List @@ e, {e}]},
    Count[t, Times[-1, __]] >
    Count[t, z_ /; !NumericQ[z] && !MatchQ[z, Times[_?NumericQ, __]]]];
DiracDelta[a___, x_, b___] /; DiracDelta`Private`flipQ[x] := DiracDelta[a, Expand[-x], b];

(* A numeric factor comes out as 1/|c|, as in the reference: DiracDelta[2 x] is
 * DiracDelta[x]/2, DiracDelta[2.5 x] is 0.4 DiracDelta[x], DiracDelta[4 x + 6 y]
 * is DiracDelta[2 x + 3 y]/2 (the rational gcd of an exact sum, constants
 * included, so DiracDelta[1 + 2 x] stays), DiracDelta[2 x, y] is
 * DiracDelta[x, y]/2.  A complex factor, or a sum with an inexact coefficient,
 * stays. *)
DiracDelta`Private`coef[t_] := Which[NumericQ[t], t, MatchQ[t, Times[_?NumericQ, __]], First[t], True, 1];
DiracDelta`Private`scale[e_] := Which[
    NumericQ[e], 1,
    Head[e] === Plus,
        With[{cs = DiracDelta`Private`coef /@ List @@ e},
            If[VectorQ[cs, MatchQ[#, _Integer | _Rational] &], GCD @@ cs, 1]],
    True,
        With[{c = DiracDelta`Private`coef[e]},
            If[NumericQ[c] && FreeQ[c, Complex] && Im[c] == 0, Abs[c], 1]]];
DiracDelta[a___, x_, b___] /; DiracDelta`Private`scale[x] =!= 1 :=
    With[{c = DiracDelta`Private`scale[x]},
        DiracDelta[a, If[Head[x] === Plus, Expand[x/c],
                         Sign[DiracDelta`Private`coef[x]] Rest[x]], b]/c];

(* Attributes as in the reference: Listable (FeynCalc's FCLoopFindSectors
 * reads sectors as HeavisideTheta[pows - 1/2] over a list), Orderless,
 * Protected.  Set last, so the rules above are installed first. *)
SetAttributes[{HeavisideTheta, DiracDelta}, {Listable, Orderless, Protected}];
