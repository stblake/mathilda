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
(* shifts integer-offset gammas. (The reflection identity Gamma[z]Gamma[1-z]   *)
(* == Pi/Sin[Pi z] is intentionally omitted: it usually raises complexity and  *)
(* emits spurious ComplexInfinity at integer z during the search.) *)
RegisterTransforms[Gamma, {
    Function[e, Replace[e, Gamma[z_ + 1] :> z Gamma[z]]]
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
