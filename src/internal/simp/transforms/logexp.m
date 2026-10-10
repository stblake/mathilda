(* ====================================================================== *)
(* logexp.m -- FullSimplify identities for logarithms of trig/hyperbolic   *)
(*             combinations that collapse to a single inverse-hyperbolic,   *)
(*             inverse-trigonometric, or linear form.                        *)
(*                                                                          *)
(* Loaded lazily by FullSimplify when the head Log appears in the input     *)
(* (every identity below has Log on the left). Each identity is a unary      *)
(* function applied (by Simplify's TransformationFunctions machinery) to     *)
(* every subexpression; it uses a top-level Replace since Simplify already   *)
(* walks every node. Identities that do not apply return their input         *)
(* unchanged.                                                                *)
(*                                                                          *)
(* THE CLASS. These are the real-log antiderivative (Gudermannian) family:   *)
(* the logarithm of a two-term +/- combination of cofunctions, or of a       *)
(* half-angle tangent/cotangent, collapses to a single inverse function.     *)
(* They have no general algorithmic route (a generic log/exp canonicaliser   *)
(* does not find them), which is exactly why they live in FullSimplify's     *)
(* rule tier rather than in Simplify. All signs were verified numerically.   *)
(*                                                                          *)
(*   Log[Sec[x] +/- Tan[x]]   == +/- ArcTanh[Sin[x]]                         *)
(*   Log[Csc[x] +/- Cot[x]]   == +/- ArcTanh[Cos[x]]                         *)
(*   Log[Cosh[x] +/- Sinh[x]] == +/- x                                       *)
(*   Log[Coth[x] +/- Csch[x]] == +/- ArcCoth[Cosh[x]]                        *)
(*   Log[Tan[x/2]]  == -ArcTanh[Cos[x]]   Log[Cot[x/2]]  == ArcTanh[Cos[x]]  *)
(*   Log[Tanh[x/2]] == -ArcCoth[Cosh[x]]  Log[Coth[x/2]] == ArcCoth[Cosh[x]] *)
(*                                                                          *)
(* TWO SURFACE FORMS PER TERM. A coefficient may reach us in either shape:   *)
(* the Simplify pipeline FACTORS a positive sum (2 Sec + 2 Tan becomes       *)
(* Times[2, Plus[Sec, Tan]]) but leaves a difference DISTRIBUTED             *)
(* (2 Sec - 2 Tan stays Plus[Times[2, Sec], Times[-2, Tan]]). So every       *)
(* two-term identity ships FOUR rules: a factored '+' and '-' with one       *)
(* optional coefficient a_. on the whole group, and a distributed '+' and    *)
(* '-' with two INDEPENDENT optional coefficients a_. (first cofunction) and *)
(* b_. (second) guarded by b === a / b === -a. The pair covers whichever     *)
(* form the pipeline hands us, and the bare case (a == 1) through either. A  *)
(* single optional coefficient on both distributed terms cannot express the  *)
(* '-' case at a non-unit coefficient, because -2 is one integer atom and    *)
(* will not match a literal -1 times c.                                      *)
(*                                                                          *)
(* POSITIVITY GUARD. Log[a z] == Log[a] + Log[z] holds on the same branch    *)
(* only for a > 0, so every rule requires TrueQ[Positive[a]]. For the bare   *)
(* forms the optional default a == 1 is positive, so they fire; a symbolic   *)
(* or negative coefficient is left untouched (sound: no branch is crossed).  *)
(* Log[a] vanishes to 0 on evaluation when a == 1. FreeQ[{a, b}, y] keeps    *)
(* the coefficients independent of the angle.                                *)
(*                                                                          *)
(* DIRECTION / TERMINATION. The rules only collapse (Log[...] -> inverse);   *)
(* no right-hand side contains a head that re-matches any left-hand side, so *)
(* there is no rewrite loop. Each collapse strictly lowers SimplifyCount, so *)
(* Simplify's search keeps it. Every statement is ';'-terminated as the file *)
(* loader requires (newlines are not statement separators).                  *)
(* ====================================================================== *)

RegisterTransforms[Log, {

    (* ------- Sec/Tan: Log[a (Sec[y] +/- Tan[y])] == Log[a] +/- ArcTanh[Sin[y]] *)
    Function[e, Replace[e,
        Log[a_. (Sec[y_] + Tan[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]]) :> Log[a] + ArcTanh[Sin[y]]]],
    Function[e, Replace[e,
        Log[a_. (Sec[y_] - Tan[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]]) :> Log[a] - ArcTanh[Sin[y]]]],
    Function[e, Replace[e,
        Log[a_. Sec[y_] + b_. Tan[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === a) :>
            Log[a] + ArcTanh[Sin[y]]]],
    Function[e, Replace[e,
        Log[a_. Sec[y_] + b_. Tan[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === -a) :>
            Log[a] - ArcTanh[Sin[y]]]],

    (* ------- Csc/Cot: Log[a (Csc[y] +/- Cot[y])] == Log[a] +/- ArcTanh[Cos[y]] *)
    Function[e, Replace[e,
        Log[a_. (Csc[y_] + Cot[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]]) :> Log[a] + ArcTanh[Cos[y]]]],
    Function[e, Replace[e,
        Log[a_. (Csc[y_] - Cot[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]]) :> Log[a] - ArcTanh[Cos[y]]]],
    Function[e, Replace[e,
        Log[a_. Csc[y_] + b_. Cot[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === a) :>
            Log[a] + ArcTanh[Cos[y]]]],
    Function[e, Replace[e,
        Log[a_. Csc[y_] + b_. Cot[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === -a) :>
            Log[a] - ArcTanh[Cos[y]]]],

    (* ------- Cosh/Sinh: Log[a (Cosh[y] +/- Sinh[y])] == Log[a] +/- y        *)
    (* (Cosh +/- Sinh == E^(+/- y); base Simplify already collapses the bare   *)
    (* case, these add the coefficient-carrying forms and insure the class.)   *)
    Function[e, Replace[e,
        Log[a_. (Cosh[y_] + Sinh[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]]) :> Log[a] + y]],
    Function[e, Replace[e,
        Log[a_. (Cosh[y_] - Sinh[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]]) :> Log[a] - y]],
    Function[e, Replace[e,
        Log[a_. Cosh[y_] + b_. Sinh[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === a) :> Log[a] + y]],
    Function[e, Replace[e,
        Log[a_. Cosh[y_] + b_. Sinh[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === -a) :> Log[a] - y]],

    (* ------- Coth/Csch: Log[a (Coth[y] +/- Csch[y])] == Log[a] +/- ArcCoth[Cosh[y]] *)
    Function[e, Replace[e,
        Log[a_. (Coth[y_] + Csch[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]]) :> Log[a] + ArcCoth[Cosh[y]]]],
    Function[e, Replace[e,
        Log[a_. (Coth[y_] - Csch[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]]) :> Log[a] - ArcCoth[Cosh[y]]]],
    Function[e, Replace[e,
        Log[a_. Coth[y_] + b_. Csch[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === a) :>
            Log[a] + ArcCoth[Cosh[y]]]],
    Function[e, Replace[e,
        Log[a_. Coth[y_] + b_. Csch[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === -a) :>
            Log[a] - ArcCoth[Cosh[y]]]],

    (* ------- Half-angle tangents/cotangents (single term; Log[k] -> 0 at k == 1). *)
    (* Log[k Tan[y/2]]  == Log[k] - ArcTanh[Cos[y]]   (real-log form of Int Csc). *)
    Function[e, Replace[e,
        Log[k_. Tan[y_/2]] /; (FreeQ[k, y] && TrueQ[Positive[k]]) :>
            Log[k] - ArcTanh[Cos[y]]]],
    (* Log[k Cot[y/2]]  == Log[k] + ArcTanh[Cos[y]].  *)
    Function[e, Replace[e,
        Log[k_. Cot[y_/2]] /; (FreeQ[k, y] && TrueQ[Positive[k]]) :>
            Log[k] + ArcTanh[Cos[y]]]],
    (* Log[k Tanh[y/2]] == Log[k] - ArcCoth[Cosh[y]]  (real-log form of Int Csch). *)
    Function[e, Replace[e,
        Log[k_. Tanh[y_/2]] /; (FreeQ[k, y] && TrueQ[Positive[k]]) :>
            Log[k] - ArcCoth[Cosh[y]]]],
    (* Log[k Coth[y/2]] == Log[k] + ArcCoth[Cosh[y]]. *)
    Function[e, Replace[e,
        Log[k_. Coth[y_/2]] /; (FreeQ[k, y] && TrueQ[Positive[k]]) :>
            Log[k] + ArcCoth[Cosh[y]]]]
}];
