(* tools/gr/check.m -- syntax self-check for one G&R transcription batch, run by
   tools/gr/verify.py --check (see tools/gr/TRANSCRIBE.md).

   grCheck[id, status, lhs, {rhs..}, {cond..}] prints one line per problem:
     @@E id message   hard error (bad status, parse failure, not exactly one
                      Integrate in lhs)
     @@W id message   warning (value names a symbol the integral lacks --
                      usually an unsubstituted abbreviation such as z, R, Delta)
   Count, not Cases, on the held lhs: Cases would evaluate the Integrate. *)

grStatuses = {"ok", "SKIP-RECURRENCE", "SKIP-MULTIPLE", "SKIP-GENERIC", "SKIP-PV",
  "SKIP-REFERENCE", "SKIP-CONTOUR", "SKIP-UNSUPPORTED", "SKIP-ILLEGIBLE"};

grCkNames[e_] := Complement[
  Union[Cases[e, grS_Symbol /; Context[grS] === "Global`" :> SymbolName[grS], {0, Infinity}, Heads -> True]],
  Union[Flatten[Cases[e, (Sum | Product | Table | Integrate)[_, grIts__List] :>
      Cases[{grIts}, {grV_Symbol, ___} :> SymbolName[grV]], {0, Infinity}, Heads -> True]]]];

grCheck[id_String, st_String, lhs_String, rhs_List, conds_List] := Module[{l, rs, cs, extra},
  If[! MemberQ[grStatuses, st], Print["@@E\t", id, "\tbad status ", st]];
  If[st === "ok",
    l = Quiet[ToExpression[lhs, InputForm, Hold]];
    If[l === $Failed, Print["@@E\t", id, "\tlhs does not parse"],
      If[Count[l, HoldPattern[Integrate[_, _]], {0, Infinity}] =!= 1,
        Print["@@E\t", id, "\tlhs must hold exactly one Integrate[f, spec]"]]];
    rs = Quiet[ToExpression[#, InputForm, Hold]] & /@ rhs;
    If[rs === {} || MemberQ[rs, $Failed], Print["@@E\t", id, "\trhs missing or does not parse"]];
    cs = Quiet[ToExpression[#, InputForm, Hold]] & /@ conds;
    If[MemberQ[cs, $Failed], Print["@@E\t", id, "\ta condition does not parse"]];
    If[l =!= $Failed && rs =!= {} && FreeQ[rs, $Failed],
      extra = Select[Complement[grCkNames[rs], grCkNames[l]], ! StringStartsQ[#, "Mathar"] &];
      If[extra =!= {}, Print["@@W\t", id, "\tvalue has symbols absent from the integral: ", extra]]]];
  Null];
