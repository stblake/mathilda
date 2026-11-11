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
(*                                                                            *)
(* Further stages (Tower, Poly, Curve, RDE, LimitedIntegrate, Hermite,        *)
(* Special, Puiseux, Residues, Divisors, Specialize, LogPart, AlgebraicPart,  *)
(* Pipeline) land phase by phase; see BRONSTEIN_PLAN.md.                      *)
(* ========================================================================== *)

LoadModule["trager/Trager.m"];
BeginPackage["Bronstein`", {"Mathilda`"}];

Begin["`Private`"];

(* Short names of Trager's private routines resolve while the siblings are    *)
(* read.  Every Bronstein definition is br-prefixed, so no Trager symbol is   *)
(* (re)defined from here.                                                     *)
$ContextPath = Append[$ContextPath, "Trager`Private`"];

LoadModule["bronstein/Common.m"];
LoadModule["bronstein/RDE.m"];

End[];
EndPackage[];
