(* ====================================================================== *)
(* FullSimplify.m -- driver and relevance engine for FullSimplify.          *)
(*                                                                          *)
(* FullSimplify[expr] tries a wide range of transformations involving       *)
(* elementary and special functions and returns the simplest form it finds. *)
(* It is a thin wrapper around Simplify: it gathers the function identities  *)
(* relevant to the heads actually present in expr and feeds them to Simplify *)
(* through the TransformationFunctions option. Because Simplify only ever    *)
(* keeps a candidate of strictly lower complexity, FullSimplify is           *)
(* guaranteed to be at least as simple as Simplify.                          *)
(*                                                                           *)
(* The full identity collection will eventually run to thousands of rules.   *)
(* To keep that scalable, identities are split into per-function-family      *)
(* libraries under simp/transforms/ that are loaded ONLY when a matching     *)
(* head appears in the input (see the manifest and Ensure below). An         *)
(* expression with no special-function heads loads no libraries and behaves  *)
(* exactly like Simplify.                                                    *)
(*                                                                           *)
(* PACKAGE LAYOUT. This is a proper Mathematica-style package: the two       *)
(* exported symbols -- FullSimplify (the user entry point) and              *)
(* RegisterTransforms (the hook the identity libraries call) -- are declared *)
(* via usage messages in the BeginPackage["FullSimplify`"] prologue; every   *)
(* helper lives in the FullSimplify`Private` context opened by               *)
(* Begin["`Private`"]. After EndPackage[], the FullSimplify` context is on   *)
(* $ContextPath, so both bare FullSimplify[...] (from any context) and the   *)
(* bare RegisterTransforms[...] calls in the lazily-loaded transform files    *)
(* resolve to the exported symbols.                                          *)
(*                                                                           *)
(* This file (the driver + manifest) is loaded eagerly at startup by         *)
(* core_init; it defines rules only and performs no simplification work at   *)
(* load time. Every top-level statement is terminated with ';' as required   *)
(* by the file loader (newlines are not statement separators).               *)
(* ====================================================================== *)

BeginPackage["FullSimplify`", {"Mathilda`"}];   (* LoadModule, SimplifyCount *)

FullSimplify::usage =
    "FullSimplify[expr] tries a wide range of transformations on expr involving elementary and special functions and returns the simplest form it finds. FullSimplify[expr, assum] simplifies using the assumptions assum. It yields at least as simple a form as Simplify. Options: ComplexityFunction, TransformationFunctions, and TimeConstraint -> {tLoc, tTot}.";

RegisterTransforms::usage =
    "RegisterTransforms[h, {f1, f2, ...}] registers unary transformation functions for head h; the FullSimplify relevance engine applies them whenever h appears in the input. Called by the identity-library files under simp/transforms/.";

Begin["`Private`"];

(* ---------------------------------------------------------------------- *)
(* Registry. Library[h] is the list of unary transformation functions       *)
(* registered for head h; unregistered heads yield {}. Each transformation   *)
(* function takes one (sub)expression and returns a candidate rewrite (or     *)
(* the input unchanged when its identity does not apply).                     *)
(* ---------------------------------------------------------------------- *)
Library[_] := {};

(* RegisterTransforms[h, {f1, f2, ...}] appends transformation functions for *)
(* head h. Called by the identity-library files as they load.                *)
RegisterTransforms[h_, fns_List] :=
    (Library[h] = Join[Library[h], fns];);

(* ---------------------------------------------------------------------- *)
(* Manifest: which identity-library file defines the transforms for a head. *)
(* Heads with no entry map to None (nothing to load). Adding a new family is *)
(* a matter of dropping a file in simp/transforms/ and listing its heads     *)
(* here. (Named ModuleFile rather than Module to avoid the System`Module     *)
(* clash.)                                                                    *)
(* ---------------------------------------------------------------------- *)
ModuleFile[_] := None;

ModuleFile[Gamma]      := "simp/transforms/gamma.m";
ModuleFile[LogGamma]   := "simp/transforms/gamma.m";
ModuleFile[PolyGamma]  := "simp/transforms/gamma.m";
ModuleFile[Pochhammer] := "simp/transforms/gamma.m";
ModuleFile[Beta]       := "simp/transforms/gamma.m";
ModuleFile[Factorial]  := "simp/transforms/gamma.m";

ModuleFile[PolyLog]    := "simp/transforms/polylog.m";

ModuleFile[Erf]        := "simp/transforms/erf.m";
ModuleFile[Erfc]       := "simp/transforms/erf.m";

ModuleFile[Surd]       := "simp/transforms/powerradical.m";

(* Log of trig/hyperbolic combinations (Gudermannian antiderivative family).
   Every identity has Log on the left, so a single Log entry loads them all. *)
ModuleFile[Log]        := "simp/transforms/logexp.m";

(* ---------------------------------------------------------------------- *)
(* Lazy library loading. Each library file is loaded at most once (guarded   *)
(* both here and, defensively, by the C loader's load-once bookkeeping).     *)
(* ---------------------------------------------------------------------- *)
FileLoaded[_] := False;

Ensure[h_] := Module[{m = ModuleFile[h]},
    If[m =!= None && FileLoaded[m] === False,
        FileLoaded[m] = True;
        LoadModule[m]
    ];
];

(* ---------------------------------------------------------------------- *)
(* Relevance engine. Collect every symbol appearing as a head or atom in     *)
(* expr, ensure each one's identity library is loaded, then return the       *)
(* de-duplicated union of all transformation functions registered for those  *)
(* heads. This is the mechanism that keeps FullSimplify from dragging in      *)
(* identities irrelevant to the input.                                       *)
(* ---------------------------------------------------------------------- *)
Relevant[expr_] := Module[{heads},
    heads = Union[Cases[expr, _Symbol, {0, Infinity}, Heads -> True]];
    Map[Ensure, heads];
    Union[Flatten[Map[Library, heads]]]
];

(* ---------------------------------------------------------------------- *)
(* Option handling. Options are trailing Rule[name, value] arguments; a      *)
(* single leading non-rule argument is the positional assumption.            *)
(* ---------------------------------------------------------------------- *)
OptValue[rules_, name_, default_] := Module[{hit},
    hit = Cases[rules, Rule[name, v_] :> v];
    If[hit === {}, default, Last[hit]]
];

(* TimeConstraint -> {tLoc, tTot} | t | Infinity, normalised to {tLoc, tTot}. *)
ParseTC[Infinity] := {Infinity, Infinity};
ParseTC[{a_, b_}] := {a, b};
ParseTC[t_]       := {t, t};

(* Wrap each transformation function with a per-transform local time budget,  *)
(* so any one transform that exceeds tLoc seconds yields its input unchanged  *)
(* (a no-op candidate) rather than stalling the whole simplification. With[]  *)
(* injects the concrete function and budget into the wrapper body: a plain     *)
(* nested Function would not close over f/tLoc.                                *)
WrapLocal[fns_, Infinity] := fns;
WrapLocal[fns_, tLoc_] :=
    Map[Function[f,
        With[{ff = f, tt = tLoc},
            Function[e, TimeConstrained[ff[e], tt, e]]]], fns];

(* ---------------------------------------------------------------------- *)
(* Separable stage. A multivariate polynomial rewritten as a minimal sum of   *)
(* separable products B[S] Q[T], each factor Factor'd                        *)
(* (simp/SeparablePolynomialDecomposition.m, loaded at startup):             *)
(*   FullSimplify[Expand[(3 a b - b^2)^3 (x^2 - y + 1)^2 - (a - b)^2 (x - y)]] *)
(*     -> -(a - b)^2 (x - y) + (3 a - b)^3 b^3 (1 + x^2 - y)^2  (LeafCount 38;  *)
(*   Simplify gives 197, Expand and Factor 251).                               *)
(* It is FullSimplify's, not Simplify's: as a Simplify stage it bought 14%    *)
(* fewer leaves on a 63-case polynomial corpus for +26% Simplify time (v0.6177; *)
(* moved here in v0.6178).                                                   *)
(*                                                                           *)
(* Which split to use cannot be known in advance, so every unordered split     *)
(* {S, T} with both sides non-empty is tried, 2^(n-1) - 1 of them, for 2 to 5  *)
(* variables (6 tripled the time for few leaves). One orientation per split:  *)
(* the last variable is always in T (both orientations: +8 leaves, +48% time). *)
(* The gate: the expanded form is an exact sum of at least 4 terms with no     *)
(* non-integer power, in symbol variables only; inexact input is left alone,  *)
(* since the decomposition of Rationalize[p, 0] is exact and FullSimplify must *)
(* not change the precision of what it is given. Only a SUM is offered: a      *)
(* rank-1 answer is a factorisation, which is Factor's (with Factor's sign     *)
(* convention, issue #126). Simplify's transformation driver keeps the        *)
(* candidate only if the whole expression gets strictly simpler.             *)
(* ---------------------------------------------------------------------- *)
SepHead = SeparablePolynomialDecomposition`SeparablePolynomialDecomposition;

SepGate[p_] := Head[p] === Plus && Length[p] >= 4 &&
    FreeQ[p, _?InexactNumberQ] && FreeQ[p, Power[_, q_ /; !IntegerQ[q]]];

SepVars[v_] := 2 <= Length[v] <= 5 && VectorQ[v, (Head[#] === Symbol)&];

(* The split numbered m (1 <= m < 2^(n-1)): S = the variables among the first
   n - 1 whose bit is set, T = the rest (always holding the last variable). *)
SepSplit[p_, vars_, m_] :=
    With[{S = Pick[Most[vars], Reverse[IntegerDigits[m, 2, Length[vars] - 1]], 1]},
        SepHead[p, S, Complement[vars, S]]];

(* SepBest[p] is memoised per FullSimplify call (Core Blocks SepMemo): the
   transformation driver re-applies every function after each improvement, and
   a node whose Expand is unchanged would redo all its splits. *)
SepBest[p_] := Module[{vars, cands},
    vars = Variables[p];
    If[!SepVars[vars], Return[None, Module]];
    cands = Select[Table[SepSplit[p, vars, m], {m, 2^(Length[vars] - 1) - 1}],
                   (Head[#] === Plus && FreeQ[#, SepHead])&];
    If[cands === {}, None, cands[[First[Ordering[Map[SimplifyCount, cands], 1]]]]]
];

Separable[e_] := Module[{p, r},
    p = Expand[e];
    If[!SepGate[p], Return[e, Module]];
    r = SepMemo[p];
    If[Head[r] === SepMemo, r = SepBest[p]; SepMemo[p] = r];
    If[r === None, e, r]
];

(* ---------------------------------------------------------------------- *)
(* Consolidation sweep. Simplify only keeps a candidate of STRICTLY lower    *)
(* complexity, so an identity that merely re-packages an expression at equal *)
(* complexity (e.g. z Gamma[z] == Gamma[z+1], where both forms have the same *)
(* SimplifyCount) is never taken by the search. This sweep gives the         *)
(* relevance-selected identities one more pass over the Simplify result and   *)
(* keeps any rewrite that does not RAISE the complexity score, consolidating  *)
(* such ties. The >= Simplify guarantee is preserved because complexity never *)
(* increases. Any strictly-simpler rewrite reachable from these identities    *)
(* was already taken during the search (they were in TransformationFunctions),*)
(* so the surviving candidates here are equicomplex re-packagings.            *)
(* ---------------------------------------------------------------------- *)
(* r, fns and score are pattern variables, so they are substituted (by value)
   throughout the body -- including inside the Select predicate and the Map
   function -- before evaluation. We therefore avoid Module locals inside those
   nested pure functions, which would NOT close over them, and recompute
   score[r] inline instead. *)
Consolidate[r_, fns_, score_] := Module[{cands},
    cands = Select[Map[Function[f, f[r]], fns],
                   (UnsameQ[#, r] && score[#] <= score[r])&];
    If[cands === {}, r, First[cands]]
];

(* ---------------------------------------------------------------------- *)
(* Driver. FullSimplify[expr, opts...] -> simplest form found.               *)
(* ---------------------------------------------------------------------- *)
FullSimplify[expr_, opts___] := Core[expr, {opts}];

Core[expr_, allOpts_List] := Module[
    {optRules, positional, assum, cf, userTF, tc, tLoc, tTot,
     relevant, allTF, simpArgs, simp, score},

    optRules   = Cases[allOpts, _Rule];
    positional = Select[allOpts, (Head[#] =!= Rule)&];

    cf     = OptValue[optRules, ComplexityFunction, Automatic];
    userTF = OptValue[optRules, TransformationFunctions, {}];
    tc     = ParseTC[OptValue[optRules, TimeConstraint, Infinity]];
    tLoc   = First[tc];
    tTot   = Last[tc];

    (* Accept either a bare function or a list of user transforms. *)
    userTF = If[Head[userTF] === List, userTF, {userTF}];

    (* Relevance-selected identities for the heads present in expr. *)
    relevant = Relevant[expr];

    (* {Automatic, Separable, <relevant>, <user-supplied>}: Automatic keeps the
       entire built-in Simplify pipeline, guaranteeing FullSimplify >= Simplify;
       Separable (above) gates itself, so it costs a non-polynomial little. *)
    allTF = Join[{Automatic}, WrapLocal[Prepend[relevant, Separable], tLoc],
                 WrapLocal[userTF, tLoc]];

    (* The positional assumption (if any) is put in force by binding
       $Assumptions = $Assumptions && assum around the whole run, exactly as
       Assuming[assum, Simplify[expr]] does -- which Simplify[expr, assum]
       is equivalent to. Binding it (rather than only forwarding it to
       Simplify) is what lets the identity libraries SEE it: a transformation
       function is an ordinary Mathilda function called during the search and
       during the consolidation sweep, and a domain-guarded identity (the
       real-log family in logexp.m) reads its facts from $Assumptions. *)
    assum = If[positional === {}, True,
        If[Head[First[positional]] === List, Apply[And, First[positional]],
           First[positional]]];

    simpArgs = Join[{expr},
                    {TransformationFunctions -> allTF, ComplexityFunction -> cf}];

    (* Final consolidation is scored with the same complexity measure Simplify
       used (SimplifyCount when ComplexityFunction is Automatic). *)
    score = If[cf === Automatic, SimplifyCount, cf];

    Block[{$Assumptions = $Assumptions && assum, SepMemo},
        (* Total time budget across all transformations. The plain-Simplify
           fallback is the HoldAll-protected third argument of TimeConstrained,
           so it is computed only on timeout -- never doubling work in the
           common (unconstrained) path -- and it preserves the >= Simplify
           contract. *)
        simp = If[tTot === Infinity,
            Apply[Simplify, simpArgs],
            TimeConstrained[
                Apply[Simplify, simpArgs],
                tTot,
                Apply[Simplify, {expr, ComplexityFunction -> cf}]]];
        Consolidate[simp, relevant, score]]
];

End[];

(* Lock the user entry point so it cannot be accidentally redefined; the
   identity-library hook RegisterTransforms stays writable so the lazily
   loaded transform files can call it. *)
SetAttributes[FullSimplify, Protected];

EndPackage[];
