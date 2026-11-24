(* latex_tex.m -- TeX for the integration report (latex_report.py).

   Run:  IN=keys.tsv OUT=tex.tsv ./Mathilda -file tools/integration_shootout/latex_tex.m

   IN has one line per expression, "key<TAB>expression" in Mathilda input
   syntax.  OUT gets "key<TAB>tex" for every line that parses, with no line for
   one that does not (the report then shows the raw text only).  RootSum[p &, g &]
   is typeset as a sum over the roots alpha of p, since TeXForm would print the
   pure functions literally; every other head is TeXForm's.  An expression is
   parsed and evaluated as written: answers are already in the form their system
   returned, and evaluating them is the same as reading them back. *)
in = Environment["IN"]; out = Environment["OUT"];
alpha = Symbol["\[Alpha]"];
texRS[e_] := Module[{rs = {}, k = 0, body, t},
  body = e /. RootSum[p_Function, g_Function] :> (k++;
     AppendTo[rs, {Symbol["rsPlaceholder" <> ToString[k]], p[alpha], g[alpha]}];
     Symbol["rsPlaceholder" <> ToString[k]]);
  t = ToString[TeXForm[body]];
  Do[t = StringReplace[t, ToString[TeXForm[r[[1]]]] ->
      "\\sum_{\\alpha:\\; " <> ToString[TeXForm[r[[2]]]] <> " = 0} " <>
      "\\left(" <> ToString[TeXForm[r[[3]]]] <> "\\right)"], {r, rs}];
  StringReplace[t, ToString[TeXForm[alpha]] -> "\\alpha"]];
lines = Select[StringSplit[Import[in, "Text"], "\n"], StringContainsQ[#, "\t"] &];
strm = OpenWrite[out];
Do[
  Module[{p = StringSplit[ln, "\t", 2], e, t},
    e = Quiet[TimeConstrained[ToExpression[p[[2]]], 20, $Failed]];
    If[e =!= $Failed && !MatchQ[e, _ToExpression],
      t = Quiet[TimeConstrained[texRS[e], 20, $Failed]];
      If[StringQ[t] && !StringContainsQ[t, "\n"], WriteString[strm, p[[1]], "\t", t, "\n"]]]],
  {ln, lines}];
Close[strm];
