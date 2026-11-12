(* ::Package:: *)

(* ========================================================================== *)
(* Bronstein.m -- Bronstein's algorithm for integrating mixed algebraic and   *)
(*                transcendental functions, restricted to ONE simple radical  *)
(*                over a Log/Exp tower.                                       *)
(*                                                                            *)
(* Integrate[f, x, Method -> "Bronstein"]   /   Integrate`Bronstein[f, x]     *)
(*                                                                            *)
(* Scope (BRONSTEIN_PLAN.md, "Scope").  Shape A: K = k(theta)(y), k = Q(x)    *)
(* (t_1..t_r) a transcendental Log/Exp tower, theta the top monomial, y^n =   *)
(* P(theta) with P depending on theta.  Shape B (radical at the bottom) is a  *)
(* later phase.  Everything else declines.                                    *)
(*                                                                            *)
(* References.                                                                *)
(*   [B]  M. Bronstein, "Integration of elementary functions", J. Symbolic    *)
(*        Computation 9 (1990) 117-173.  Page numbers are the JOURNAL ones    *)
(*        (the PDF page is the journal page - 116).                           *)
(*   [T]  B. M. Trager, "Integration of Algebraic Functions", Ph.D. thesis,   *)
(*        MIT, 1984, Appendix A (printed page numbers).                       *)
(*                                                                            *)
(* Architecture: one file per stage, every function br-prefixed (no clash     *)
(* with Trager's names) and carrying a ::usage docstring.  The Trager package *)
(* is loaded first and Trager`Private` is put on $ContextPath while the       *)
(* siblings are read, so its field / ideal / torsion routines (kc, hnf,       *)
(* principalGenerator, ...) are called by their short names.  No file assigns *)
(* a top-level variable (it would land in Global`): configuration values are  *)
(* zero-argument functions (brOption).                                        *)
(*                                                                            *)
(*   Common.m      failure protocol (Trager's), options                       *)
(*   RDE.m         tri-state Risch DE over the tower (Risch`RischDECertified) *)
(*   Tower.m       surface kernels -> monomials, radical, shape, flattening   *)
(*   Poly.m        k[theta]: content, normal/special split, derivation        *)
(*   Curve.m       canonical form over k(theta), component curves Y^m = G     *)
(*   LimitedIntegrate.m  integration in k: parametric / limited integration *)
(*   Hermite.m     Hermite reduction at the normal places                     *)
(*   AlgebraicPart.m  D b + s b = c over k(theta) (certified)                 *)
(*   Puiseux.m     expansions at theta = 0 and infinity (binomial series)     *)
(*   Special.m     reduction at the special places (exponential theta)        *)
(*   Residues.m    R(z) and the resultant criterion; special-place residues   *)
(*   Divisors.m    residue divisors over k(t), principal generators           *)
(*   Specialize.m  torsion order by specialization of k                       *)
(*   LogPart.m     log part, weight-1 projection, remainder                   *)
(*   Pipeline.m    components, reassembly, verification, BronsteinIntegrate  *)
(*                                                                            *)
(*   K0.m          Shape B: arithmetic in k0 = k(Y) and k0[theta]             *)
(*   BaseAlgebraic.m  Shape B, the component method (radical below the tower) *)
(* ========================================================================== *)

LoadModule["trager/Trager.m"];
BeginPackage["Bronstein`", {"Mathilda`"}];

BronsteinIntegrate::usage = "BronsteinIntegrate[f, x] integrates f, an elementary function of x built from rational operations, logarithms, exponentials and ONE radical (n-th root) whose radicand involves the top transcendental monomial (Bronstein 1990). It returns an antiderivative (verified), {\"not elementary\", reason} (a certificate), or {\"failed\", reason} (out of scope or undecided).";

Begin["`Private`"];

(* Short names of Trager's public and private routines resolve while the     *)
(* siblings are read.  Every Bronstein definition is br-prefixed, so no Trager *)
(* symbol is (re)defined from here.  A public name of another package that is *)
(* NOT on this path would silently become a new, undefined Bronstein`Private` *)
(* symbol (an unevaluated call): qualify such names (Risch`..., Integrate`...) *)
(* -- the suites assert that no capitalised symbol lands in this context.    *)
$ContextPath = Join[$ContextPath, {"Trager`", "Trager`Private`"}];

LoadModule["bronstein/Common.m"];
LoadModule["bronstein/RDE.m"];
LoadModule["bronstein/Tower.m"];
LoadModule["bronstein/Poly.m"];
LoadModule["bronstein/Curve.m"];
LoadModule["bronstein/LimitedIntegrate.m"];
LoadModule["bronstein/Hermite.m"];
LoadModule["bronstein/AlgebraicPart.m"];
LoadModule["bronstein/Puiseux.m"];
LoadModule["bronstein/Special.m"];
LoadModule["bronstein/Residues.m"];
LoadModule["bronstein/Divisors.m"];
LoadModule["bronstein/Specialize.m"];
LoadModule["bronstein/LogPart.m"];
LoadModule["bronstein/K0.m"];
LoadModule["bronstein/BaseAlgebraic.m"];
LoadModule["bronstein/Pipeline.m"];

End[];
EndPackage[];
