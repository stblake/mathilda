(* ::Package:: *)

(* ========================================================================== *)
(* Trager.m -- Trager's algorithm for integrating algebraic functions over a  *)
(*             simple radical extension, with Kauers' heuristic as an          *)
(*             alternative engine for the logarithmic part.                   *)
(*                                                                            *)
(* Integrate[f, x, Method -> "Trager"]        (complete algorithm)            *)
(* Integrate[f, x, Method -> "TragerKauers"]  (Kauers 2008 log-part heuristic)*)
(*                                                                            *)
(* Scope.  f is a rational function of x and ONE radical r(x)^(1/n), r in     *)
(* Q(x), with rational coefficients (several powers of the same radicand are  *)
(* allowed: Sqrt[x] and x^(1/3) are both powers of x^(1/6)).                  *)
(*                                                                            *)
(* References.                                                                *)
(*   [T]  B. M. Trager, "Integration of Algebraic Functions", Ph.D. thesis,   *)
(*        MIT, 1984.  Section and page numbers below are the PRINTED ones.    *)
(*   [B]  M. Bronstein, "Symbolic Integration Tutorial", ISSAC'98.            *)
(*   [K]  M. Kauers, "Integration of algebraic functions: a simple heuristic  *)
(*        for finding the logarithmic part", ISSAC 2008.                      *)
(*                                                                            *)
(* Architecture (one file per stage; every function carries a ::usage         *)
(* docstring).                                                                *)
(*                                                                            *)
(*   Common.m      failure protocol, small predicates                         *)
(*   Field.m       coefficient fields Q, Q(theta), GF(p) and polynomial       *)
(*                 arithmetic over them (division, extended gcd, linear       *)
(*                 algebra)                                                   *)
(*   Validate.m    radical detection, the canonical form f = Sum a_i(x) y^i   *)
(*   Reduce.m      absolute irreducibility of Y^m - G     [T] Ch 3 Sec 3 p.39 *)
(*   Genus.m       Riemann-Hurwitz genus                  [T] Ch 2 Sec 4 p.29 *)
(*   Basis.m       integral basis w_i = Y^i/d_i           [T] Ch 2 Sec 5 p.30 *)
(*   Arithmetic.m  arithmetic in K(x)[Y]/(Y^m - G) on that basis              *)
(*   Shift.m       move infinity to a finite point        [T] Ch 4 Sec 1 p.44 *)
(*   Hermite.m     algebraic Hermite reduction            [T] Ch 4 Sec 2 p.44 *)
(*   Residues.m    residues, their Q-basis, residue field [T] Ch 5 Sec 2 p.56 *)
(*   Normalize.m   Hermite normal form over K[x] and normalisation at         *)
(*                 infinity                               [T] Ch 2 Sec 3 p.27 *)
(*   Ideal.m       fractional ideals = divisors           [T] Ch 5 Sec 3 p.59 *)
(*   PrincipalGen.m principality test                     [T] Ch 6 Sec 1 p.64 *)
(*   TorsionModP.m order of a divisor class by reduction  [T] Ch 6 Sec 2-3    *)
(*                 modulo two primes                          pp.66-70        *)
(*   TragerLogTerms.m  logarithmic part, Trager's construction                *)
(*   KauersLogTerms.m  logarithmic part, Kauers' Groebner heuristic           *)
(*   LogToArcTan.m conjugate logarithm pairs -> ArcTan                        *)
(*   Pipeline.m    the per-component driver and the public entry point        *)
(*                                                                            *)
(* Design: the y^i components are integrated separately.  For a simple       *)
(* radical, f = Sum_i a_i(x) y^i is elementary iff every a_i y^i dx is, and   *)
(* the algebraic part of the integral of a_i y^i is v_i(x) y^i ([T] Appendix  *)
(* A, p.74).  Each component a_i y^i is a rational function times the single  *)
(* generator Y = y^i/s(x) of the subfield Q(x, y^i), with Y^m = G(x),         *)
(* G squarefree-ish (every exponent < m) and absolutely irreducible.  Working  *)
(* there, every pole of the Hermite remainder at a branch place has residue   *)
(* zero ([T] Ch 5 Sec 4, pp.62-63), so every residue divisor is supported on  *)
(* unramified places and is the ideal (N Y - alpha E', B) of the integral     *)
(* closure, exactly.                                                          *)
(* ========================================================================== *)

BeginPackage["Trager`", {"Mathilda`"}];

TragerIntegrate::usage = "TragerIntegrate[f, x, method] integrates f, a rational function of x and powers of a single radical r(x)^(1/n) with rational coefficients, by Trager's algorithm for simple radical extensions. method is \"Trager\" (the complete algorithm: divisors, principal generators, torsion orders by reduction modulo primes) or \"Kauers\" (Kauers' Groebner-basis heuristic for the logarithmic part). It returns
  an antiderivative                    -- verified exactly before it is returned;
  {\"not elementary\", reason}         -- a certificate that the integral is not elementary;
  {\"failed\", reason}                 -- no verdict (out of scope, or an incomplete heuristic).";

TragerIntegrateCurve::usage = "TragerIntegrateCurve[a, x, m, G, method] integrates the single component a(x) Y dx on the curve Y^m = G(x), where G is a polynomial over Q whose irreducible factors all have multiplicity < m and Y^m - G is absolutely irreducible. It returns <|\"Algebraic\" -> v, \"Logs\" -> {{c, u}, ...}|> with the integral v + Sum c Log[u] (in x and the symbol Y), or a list {\"not elementary\", reason} / {\"failed\", reason} as TragerIntegrate does.";

Begin["`Private`"];

(* Siblings are read INTO this private context, so their short function names *)
(* resolve to each other (a short private function name does not resolve      *)
(* through $ContextPath).  Order matters only for top-level side effects;      *)
(* none of these files has any.  No file assigns a top-level variable: a       *)
(* top-level variable in a LoadModule'd file lands in Global` (functions do    *)
(* not), so configuration values are zero-argument functions instead.          *)
LoadModule["trager/Common.m"];
LoadModule["trager/Field.m"];
LoadModule["trager/Validate.m"];
LoadModule["trager/Reduce.m"];
LoadModule["trager/Genus.m"];
LoadModule["trager/Basis.m"];
LoadModule["trager/Arithmetic.m"];
LoadModule["trager/Shift.m"];
LoadModule["trager/Hermite.m"];
LoadModule["trager/Residues.m"];
LoadModule["trager/Normalize.m"];
LoadModule["trager/Ideal.m"];
LoadModule["trager/PrincipalGen.m"];
LoadModule["trager/TorsionModP.m"];
LoadModule["trager/TragerLogTerms.m"];
LoadModule["trager/KauersLogTerms.m"];
LoadModule["trager/LogToArcTan.m"];
LoadModule["trager/Pipeline.m"];

End[];
EndPackage[];
