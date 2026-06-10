(* ====================================================================== *)
(* powerradical.m -- FullSimplify identities for real radicals.            *)
(*                                                                          *)
(* Loaded lazily when Surd appears in the input. Surd is the real-valued    *)
(* n-th root, so Surd[x, n]^n == x holds for every real x -- a reduction the *)
(* evaluator and Simplify do not perform on their own (unlike Sqrt[x]^2,     *)
(* which auto-evaluates to x and therefore needs no rule here).              *)
(* ====================================================================== *)

(* Surd[x, n]^n == x. The exponent pattern reuses the index variable n_, so  *)
(* it matches only when the power equals the root index. *)
RegisterTransforms[Surd, {
    Function[e, Replace[e, Surd[x_, n_]^n_ :> x]]
}];
