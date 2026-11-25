(* normalize_constants.m -- remove the constant factors from a corpus's integrands
   (stripConstants, common.m), for normalize_constants.py.

   Run:  IN=in.tsv OUT=out.tsv ./Mathilda -file tools/integration_shootout/normalize_constants.m

   IN: "CASE<TAB>id<TAB>class<TAB>family<TAB>expect<TAB>f<TAB>F<TAB>points<TAB>values<TAB>fullform",
   the corpus line as it stands.  OUT gets one line per case: "SAME<TAB>id" when f has no
   constant factor, else "NEW<TAB>id<TAB>verdict<TAB>" and the rewritten corpus line, with
   f/k for f, F/k for the reference (F is None stays None) and the values at the same
   points recomputed.  verdict is "ok" when k is free of x and, for a known F,
   PossibleZeroQ[D[F/k, x] - f/k] is True; anything else is reported and the case kept
   as it was. *)
Get["tools/integration_shootout/common.m"];
in = Environment["IN"]; out = Environment["OUT"];
put[s_String] := Module[{strm = OpenAppend[out]}, WriteString[strm, s, "\n"]; Close[strm]];
Do[
  Module[{p = StringSplit[ln, "\t"], f, F, g, k, Fn, pts, fv, verdict},
    f = ToExpression[p[[6]]]; F = If[p[[7]] === "None", None, ToExpression[p[[7]]]];
    pts = ToExpression[p[[8]]];
    {k, g} = TimeConstrained[stripConstantsK[f], 30, {1, f}];
    If[g === f, put["SAME\t" <> p[[2]]],
      verdict = Which[
        !FreeQ[k, x] || !NumericQ[k], "k not constant",
        !TrueQ[TimeConstrained[PossibleZeroQ[f - k g], 60, False]], "f != k g",
        F === None, Fn = None; "ok",
        True,
          Fn = scDivide[F, k];
          Replace[TimeConstrained[PossibleZeroQ[D[Fn, x] - g], 60, $Aborted],
            {True -> "ok", False -> "PossibleZeroQ False", _ -> "PossibleZeroQ undecided"}]];
      If[verdict =!= "ok", put["NEW\t" <> p[[2]] <> "\t" <> verdict <> "\t" <> ln],
        fv = N[g /. x -> #, 40] & /@ pts;
        put[StringRiffle[{"NEW", p[[2]], "ok", "CASE", p[[2]], p[[3]], p[[4]], p[[5]],
          ToString[InputForm[g]], If[Fn === None, "None", ToString[InputForm[Fn]]],
          p[[8]], ToString[InputForm[fv]], ToString[FullForm[g]]}, "\t"]]]]],
  {ln, Select[StringSplit[Import[in, "Text"], "\n"], StringStartsQ[#, "CASE\t"] &]}];
