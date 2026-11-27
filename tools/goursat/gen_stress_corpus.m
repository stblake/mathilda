(* tools/goursat/gen_stress_corpus.m -- corpus.tsv of benchmark 265-goursat-stress.

     ./Mathilda -file tools/goursat/gen_stress_corpus.m > benchmarks/265-goursat-stress/corpus.tsv

   One CASE row per case of every GeneralisedGoursat corpus in
   tests/data/generalised_goursat/ (unit, examples, literature, welz, stress_2026,
   stress_7, stress_11, pseudoelliptic_sqrt, pseudoelliptic_cbrt), in the shootout's
   format (tools/integration_shootout/shootout.py load_corpus):

     CASE id class family expect f F points f(points) fullform | corpus label truth

   class is "G", family the corpus's own Family field (sqrt-klein, cbrt-mixed,
   pseudoelliptic-sqrt, ...), expect elem / nonelem / unknown from the corpus's Truth
   ("elementary" / "non-elementary" / "undecided" or None), F None (the corpora carry
   no reference antiderivative).  The integrand is the corpus's own, with t renamed
   x, and is NOT stripped of constant factors: it must stay the integrand goursat.py
   was timed on (tests/data/generalised_goursat/sympy_*.tsv).  The three columns after
   fullform -- corpus, label, truth -- tie a row back to the corpus and to those
   timings.

   Points: three of the shootout's rational candidates (common.m) where f is finite
   and, preferably, real with positive radicands; failing that, any where f is a
   finite number (some integrands have complex coefficients).  They feed only the
   systems' own point checks, which decide nothing. *)

ggCorpora = {"unit", "examples", "literature", "welz", "stress_2026", "stress_7",
             "stress_11", "pseudoelliptic_sqrt", "pseudoelliptic_cbrt"};
ggShort = <|"unit" -> "unit", "examples" -> "ex", "literature" -> "lit", "welz" -> "welz",
            "stress_2026" -> "s2026", "stress_7" -> "s7", "stress_11" -> "s11",
            "pseudoelliptic_sqrt" -> "psqrt", "pseudoelliptic_cbrt" -> "pcbrt"|>;
cands = {3/7, 5/4, 7/3, 13/5, 2/9, 11/3, 17/4, 6/5, 29/7, 1/11, 23/3, 31/9, 9/7, 8/3};

finiteAt[f_, p_] := Module[{v = Quiet[N[f /. x -> p, 30]]}, NumberQ[v] && Abs[v] < 10^8];
realAt[f_, p_] := Module[{bases, v},
  bases = Cases[{f}, Power[b_, _Rational] :> b, {0, Infinity}];
  v = Quiet[N[f /. x -> p, 30]];
  NumberQ[v] && Im[v] == 0 && Abs[v] < 10^8 &&
    AllTrue[Quiet[N[bases /. x -> p, 30]], NumberQ[#] && Im[#] == 0 && # > 1/100 &]];
ggPoints[f_] := Module[{good = Select[cands, realAt[f, #] &]},
  If[Length[good] < 3, good = Join[good, Select[Complement[cands, good], finiteAt[f, #] &]]];
  Take[good, UpTo[3]]];

ggExpect[tr_] := Switch[tr, "elementary", "elem", "non-elementary", "nonelem", _, "unknown"];
ggTruth[tr_] := If[StringQ[tr], tr, "None"];
ggClean[s_String] := StringReplace[s, {"\t" -> " ", "\n" -> " "}];

Do[
  Module[{cases = Get["tests/data/generalised_goursat/" <> corp <> ".m"], k = <||>},
    Do[
      Module[{f = c["Integrand"] /. c["Var"] -> x, id, pts, fv},
        (* numbered within its family, so G-s7-sqrt-klein-004 is "sqrt-klein #4" *)
        id = "G-" <> ggShort[corp] <> "-" <> c["Family"] <> "-" <>
          IntegerString[Lookup[k, c["Family"], 0], 10, 3];
        k[c["Family"]] = Lookup[k, c["Family"], 0] + 1;
        pts = TimeConstrained[ggPoints[f], 30, {}];
        fv = N[f /. x -> #, 40] & /@ pts;
        Print["CASE\t", id, "\tG\t", c["Family"], "\t", ggExpect[c["Truth"]], "\t",
          ToString[InputForm[f]], "\tNone\t", ToString[InputForm[pts]], "\t",
          ToString[InputForm[fv]], "\t", ToString[FullForm[f]], "\t", corp, "\t",
          ggClean[c["Label"]], "\t", ggTruth[c["Truth"]]]],
      {c, cases}]],
  {corp, ggCorpora}];
