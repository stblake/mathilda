(* ====================================================================== *)
(* erf.m -- FullSimplify identities for the error-function family.          *)
(*                                                                          *)
(* Loaded lazily when Erf, Erfc, or Erfi appears in the input. Parity        *)
(* (Erf[-z] -> -Erf[z], etc.) is already handled by the kernel's builtins,   *)
(* so the value FullSimplify adds here is the complementary identity that    *)
(* Simplify does not perform on its own.                                     *)
(* ====================================================================== *)

(* Erf[z] + Erfc[z] == 1  (complementary error function). The Orderless Plus *)
(* match also collapses the pair inside a larger sum, e.g. a + Erf[x] +       *)
(* Erfc[x] -> 1 + a.                                                          *)
RegisterTransforms[Erf, {
    Function[e, Replace[e, Erf[z_] + Erfc[z_] :> 1]]
}];

(* Erfc, expressed via Erf, so the complementary identity is reachable        *)
(* whichever head the expression is collected under. *)
RegisterTransforms[Erfc, {
    Function[e, Replace[e, Erf[z_] + Erfc[z_] :> 1]]
}];
