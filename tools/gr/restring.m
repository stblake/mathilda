(* restring.m -- re-print the strings of an oracle record (no integration).

   Driven by tools/gr/restring.py.  For each record it re-runs the oracle's
   own isolation (verify.m grIsolate) under the record's elliptic placeholder
   variant and prints the integrand, variable, limits and value of every rhs
   alternative with the CURRENT printer:

     @@S  id  var  integrand  limits  rhsIndex  value

   The oracle wrote these strings with ToString[..., InputForm], and before
   v0.6488 that printed a held -1/q^(n + 1) as -q^(-n + 1): the verdicts were
   sound (the oracle evaluates the held expression), but the stored strings
   that feed the corpus and the GR table could be a different expression. *)

Get["tools/gr/verify.m"];

(* Every name here is gr-prefixed, as in verify.m: the integrand is
   evaluated with the book's own symbols, and a Do iterator named k (or a
   Module local the book also uses) would substitute into it. *)
grRestring[id_String, lhs_String, rhs_List, variant_Integer] := Module[
  {lhsH, rhsL, vr, iso, f, spec, var, lim},
  lhsH = Quiet[ToExpression[lhs, InputForm, Hold]];
  rhsL = Quiet[ToExpression[#, InputForm, Hold]] & /@ rhs;
  vr = If[grHasPlaceholder[{lhsH, rhsL}] && 1 <= variant <= Length[grVariants],
          grVariants[[variant]], {}];
  Do[
    If[rhsL[[grRk]] === $Failed, Continue[]];
    iso = grIsolate[lhsH /. vr, rhsL[[grRk]] /. vr];
    If[iso === $Failed, Continue[]];
    f = Extract[iso[[1]], 1]; spec = Extract[iso[[1]], 2];
    var = If[ListQ[spec], spec[[1]], spec];
    lim = If[ListQ[spec], Rest[spec], {}];
    Print["@@S\t", id, "\t", grStr[var], "\t", grStr[f], "\t",
      If[lim === {}, "", grStr[lim]], "\t", grRk, "\t", grStrHeld[iso[[2]]]],
    {grRk, Length[rhsL]}];
  Null];
