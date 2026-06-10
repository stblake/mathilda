(* ====================================================================== *)
(* FullSimplify.m -- driver and relevance engine for FullSimplify.         *)
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
(* head appears in the input (see the manifest and FullSimplify`Ensure       *)
(* below). An expression with no special-function heads loads no libraries   *)
(* and behaves exactly like Simplify.                                        *)
(*                                                                           *)
(* This file (the driver + manifest) is loaded eagerly at startup by         *)
(* core_init; it defines rules only and performs no simplification work at   *)
(* load time. Every top-level statement is terminated with ';' as required   *)
(* by the file loader (newlines are not statement separators).               *)
(* ====================================================================== *)

(* ---------------------------------------------------------------------- *)
(* Registry. FullSimplify`Library[h] is the list of unary transformation   *)
(* functions registered for head h; unregistered heads yield {}. Each       *)
(* transformation function takes one (sub)expression and returns a candidate *)
(* rewrite (or the input unchanged when its identity does not apply).        *)
(* ---------------------------------------------------------------------- *)
FullSimplify`Library[_] := {};

(* RegisterTransforms[h, {f1, f2, ...}] appends transformation functions for *)
(* head h. Called by the identity-library files as they load.                *)
RegisterTransforms[h_, fns_List] :=
    (FullSimplify`Library[h] = Join[FullSimplify`Library[h], fns];);

(* ---------------------------------------------------------------------- *)
(* Manifest: which identity-library file defines the transforms for a head. *)
(* Heads with no entry map to None (nothing to load). Adding a new family is *)
(* a matter of dropping a file in simp/transforms/ and listing its heads     *)
(* here.                                                                      *)
(* ---------------------------------------------------------------------- *)
FullSimplify`Module[_] := None;

FullSimplify`Module[Gamma]      := "simp/transforms/gamma.m";
FullSimplify`Module[LogGamma]   := "simp/transforms/gamma.m";
FullSimplify`Module[PolyGamma]  := "simp/transforms/gamma.m";
FullSimplify`Module[Pochhammer] := "simp/transforms/gamma.m";
FullSimplify`Module[Beta]       := "simp/transforms/gamma.m";
FullSimplify`Module[Factorial]  := "simp/transforms/gamma.m";

FullSimplify`Module[PolyLog]    := "simp/transforms/polylog.m";

FullSimplify`Module[Erf]        := "simp/transforms/erf.m";
FullSimplify`Module[Erfc]       := "simp/transforms/erf.m";

FullSimplify`Module[Surd]       := "simp/transforms/powerradical.m";

(* ---------------------------------------------------------------------- *)
(* Lazy library loading. Each library file is loaded at most once (guarded   *)
(* both here and, defensively, by the C loader's load-once bookkeeping).     *)
(* ---------------------------------------------------------------------- *)
FullSimplify`FileLoaded[_] := False;

FullSimplify`Ensure[h_] := Module[{m = FullSimplify`Module[h]},
    If[m =!= None && FullSimplify`FileLoaded[m] === False,
        FullSimplify`FileLoaded[m] = True;
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
FullSimplify`Relevant[expr_] := Module[{heads},
    heads = Union[Cases[expr, _Symbol, {0, Infinity}, Heads -> True]];
    Map[FullSimplify`Ensure, heads];
    Union[Flatten[Map[FullSimplify`Library, heads]]]
];

(* ---------------------------------------------------------------------- *)
(* Option handling. Options are trailing Rule[name, value] arguments; a      *)
(* single leading non-rule argument is the positional assumption.            *)
(* ---------------------------------------------------------------------- *)
FullSimplify`OptValue[rules_, name_, default_] := Module[{hit},
    hit = Cases[rules, Rule[name, v_] :> v];
    If[hit === {}, default, Last[hit]]
];

(* TimeConstraint -> {tLoc, tTot} | t | Infinity, normalised to {tLoc, tTot}. *)
FullSimplify`ParseTC[Infinity] := {Infinity, Infinity};
FullSimplify`ParseTC[{a_, b_}] := {a, b};
FullSimplify`ParseTC[t_]       := {t, t};

(* Wrap each transformation function with a per-transform local time budget,  *)
(* so any one transform that exceeds tLoc seconds yields its input unchanged  *)
(* (a no-op candidate) rather than stalling the whole simplification. With[]  *)
(* injects the concrete function and budget into the wrapper body: a plain     *)
(* nested Function would not close over f/tLoc.                                *)
FullSimplify`WrapLocal[fns_, Infinity] := fns;
FullSimplify`WrapLocal[fns_, tLoc_] :=
    Map[Function[f,
        With[{ff = f, tt = tLoc},
            Function[e, TimeConstrained[ff[e], tt, e]]]], fns];

(* ---------------------------------------------------------------------- *)
(* Driver. FullSimplify[expr, opts...] -> simplest form found.               *)
(* ---------------------------------------------------------------------- *)
FullSimplify[expr_, opts___] := FullSimplify`Core[expr, {opts}];

FullSimplify`Core[expr_, allOpts_List] := Module[
    {optRules, positional, assumSeq, cf, userTF, tc, tLoc, tTot,
     relevant, allTF, simpArgs},

    optRules   = Cases[allOpts, _Rule];
    positional = Select[allOpts, (Head[#] =!= Rule)&];

    cf     = FullSimplify`OptValue[optRules, ComplexityFunction, Automatic];
    userTF = FullSimplify`OptValue[optRules, TransformationFunctions, {}];
    tc     = FullSimplify`ParseTC[
                 FullSimplify`OptValue[optRules, TimeConstraint, Infinity]];
    tLoc   = First[tc];
    tTot   = Last[tc];

    (* Accept either a bare function or a list of user transforms. *)
    userTF = If[Head[userTF] === List, userTF, {userTF}];

    (* Relevance-selected identities for the heads present in expr. *)
    relevant = FullSimplify`Relevant[expr];

    (* {Automatic, <relevant>, <user-supplied>}: Automatic keeps the entire
       built-in Simplify pipeline, guaranteeing FullSimplify >= Simplify. *)
    allTF = Join[{Automatic}, FullSimplify`WrapLocal[relevant, tLoc],
                 FullSimplify`WrapLocal[userTF, tLoc]];

    (* Forward the positional assumption (if any) so Simplify combines it with
       $Assumptions exactly as it does for Simplify[expr, assum]. *)
    assumSeq = If[positional === {}, {}, {First[positional]}];

    simpArgs = Join[{expr}, assumSeq,
                    {TransformationFunctions -> allTF, ComplexityFunction -> cf}];

    (* Total time budget across all transformations. The plain-Simplify
       fallback is the HoldAll-protected third argument of TimeConstrained, so
       it is computed only on timeout -- never doubling work in the common
       (unconstrained) path -- and it preserves the >= Simplify contract. *)
    If[tTot === Infinity,
        Apply[Simplify, simpArgs],
        TimeConstrained[
            Apply[Simplify, simpArgs],
            tTot,
            Apply[Simplify, Join[{expr}, assumSeq, {ComplexityFunction -> cf}]]]]
];
