(* ====================================================================== *)
(* gamma.m -- FullSimplify identities for the gamma-function family.       *)
(*                                                                          *)
(* Loaded lazily by FullSimplify when one of the heads Gamma, LogGamma,     *)
(* PolyGamma, Pochhammer, Beta, or Factorial appears in the input. Each      *)
(* identity is a unary function applied (by Simplify's TransformationFunc-   *)
(* tions machinery) to every subexpression; it uses top-level Replace since  *)
(* Simplify already walks every node. Identities that do not apply return    *)
(* their input unchanged.                                                    *)
(* ====================================================================== *)

(* Gamma[z+1] == z Gamma[z]  (recurrence; the engine reuses it both to       *)
(* collapse Gamma[x+1]/Gamma[x] -> x and to lower integer-shifted gammas).   *)
(* Recurrence Gamma[z+1] == z Gamma[z]: collapses Gamma[x+1]/Gamma[x] -> x and *)
(* shifts integer-offset gammas.                                              *)
(*                                                                            *)
(* The reflection / conjugate-pair identities below fire ONLY as PAIR rules   *)
(* (both gammas must be present in the same Times), so they never raise the   *)
(* gamma count, and each is guarded to avoid the spurious ComplexInfinity at  *)
(* integer z that made a naive lone-gamma reflection unsafe during the search.*)
RegisterTransforms[Gamma, {
    Function[e, Replace[e, Gamma[z_ + 1] :> z Gamma[z]]],
    (* Raising / consolidation: z Gamma[z] == Gamma[z+1]. Absorbs a factor   *)
    (* that matches the gamma argument (e.g. x Gamma[x] -> Gamma[1+x]). The   *)
    (* Orderless Times match lets the factor and the Gamma appear in any      *)
    (* position, and the leading r___ carries any unrelated cofactors along.  *)
    Function[e, Replace[e,
        HoldPattern[Times[r___, z_, Gamma[z_]]] :> Gamma[z + 1] Times[r]]],
    (* Reflection PAIR: Gamma[z] Gamma[1-z] == Pi / Sin[Pi z].  Fires only    *)
    (* when the two arguments sum to 1 and z is not an integer (guard against *)
    (* Sin[Pi n] == 0 -> ComplexInfinity).  e.g. Gamma[1/4] Gamma[3/4] ->     *)
    (* Pi Sqrt[2].                                                            *)
    Function[e, Replace[e,
        HoldPattern[Times[r___, Gamma[z_], Gamma[w_]]] /;
            (! IntegerQ[z]) && PossibleZeroQ[z + w - 1] :>
            Pi / Sin[Pi z] Times[r]]],
    (* Conjugate PAIR (Re == 1): Gamma[1+I b] Gamma[1-I b] == Pi b/Sinh[Pi b].*)
    Function[e, Replace[e,
        HoldPattern[Times[r___, Gamma[u_], Gamma[v_]]] /;
            (! PossibleZeroQ[Im[u]]) && PossibleZeroQ[u - Conjugate[v]]
            && PossibleZeroQ[Re[u] - 1] :>
            (Pi Im[u] / Sinh[Pi Im[u]]) Times[r]]],
    (* Conjugate PAIR (Re == 1/2): Gamma[1/2+I b] Gamma[1/2-I b] == Pi/Cosh[Pi b].*)
    Function[e, Replace[e,
        HoldPattern[Times[r___, Gamma[u_], Gamma[v_]]] /;
            (! PossibleZeroQ[Im[u]]) && PossibleZeroQ[u - Conjugate[v]]
            && PossibleZeroQ[Re[u] - 1/2] :>
            (Pi / Cosh[Pi Im[u]]) Times[r]]]
}];

(* LogGamma[z+1] == Log[z] + LogGamma[z]. *)
RegisterTransforms[LogGamma, {
    Function[e, Replace[e, LogGamma[z_ + 1] :> Log[z] + LogGamma[z]]]
}];

(* PolyGamma[n, z+1] == PolyGamma[n, z] + (-1)^n n! / z^(n+1). *)
RegisterTransforms[PolyGamma, {
    Function[e, Replace[e,
        PolyGamma[n_, z_ + 1] :> PolyGamma[n, z] + (-1)^n n! / z^(n + 1)]],
    (* PolyGamma[0, z+1] == PolyGamma[0, z] + 1/z. *)
    Function[e, Replace[e,
        PolyGamma[0, z_ + 1] :> PolyGamma[0, z] + 1/z]]
}];

(* Pochhammer[a, n] == Gamma[a+n]/Gamma[a]. *)
RegisterTransforms[Pochhammer, {
    Function[e, Replace[e, Pochhammer[a_, n_] :> Gamma[a + n]/Gamma[a]]]
}];

(* Beta[a, b] == Gamma[a] Gamma[b] / Gamma[a+b]. *)
RegisterTransforms[Beta, {
    Function[e, Replace[e, Beta[a_, b_] :> Gamma[a] Gamma[b] / Gamma[a + b]]]
}];

(* Factorial[n] == Gamma[n+1]; lets gamma identities reach factorials. *)
RegisterTransforms[Factorial, {
    Function[e, Replace[e, Factorial[n_] :> Gamma[n + 1]]]
}];
