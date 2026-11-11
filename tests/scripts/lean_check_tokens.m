(* lean_check_tokens.m -- run by lean_check_tokens.sh; see there.

   Reads the output of tools/lean_command_tokens.lean (MATHILDA_TOKENS_FILE)
   and checks, against Lean`Internal`CheckSource itself (so the property
   tested is the checker's behaviour, not a copy of its list):

   1. every leading token of every command parser is refused when it occurs
      inside a command (a continuation line), and also at column 0 unless it
      is one of the allowed command starters;
   2. no leading command parser is unindexed (it could start at an
      identifier, which no word list could catch);
   3. for every token of the token table that ends in a word character
      [A-Za-z0-9_] but is not identifier-shaped (Mathlib's `->o`, `[MOD`, ...;
      Lean's next token starts right after it, inside one of the checker's
      words), a keyword glued to it is still refused.

   MATHILDA_TOKENS_MATHLIB=1 checks with mathlib -> True.  MATHILDA_TOKENS_WORDS=1
   prints, instead, the [A-Za-z0-9_] words of every leading token that is not
   refused lexically: the contents of src/lean_check_tokens.h. *)

file = Environment["MATHILDA_TOKENS_FILE"];
ml = Environment["MATHILDA_TOKENS_MATHLIB"] === "1";
lines = StringSplit[ReadString[file], "\n"];
lead = DeleteDuplicates[StringDrop[#, 5] & /@ Select[lines, StringStartsQ[#, "lead "] &]];
toks = DeleteDuplicates[StringDrop[#, 6] & /@ Select[lines, StringStartsQ[#, "token "] &]];
nonidx = Select[lines, StringStartsQ[#, "nonindexed "] &];

(* Allowed command starters (column 0 only). *)
starters = {"theorem", "lemma", "example", "open", "namespace", "section",
  "noncomputable", "end", "variable", "universe", "set_option", "import"};
(* Leading tokens that cannot start a command on their own, with the reason. *)
exempt = <|"?" -> "only before binder_predicate (a refused word) or after an antiquotation $ (refused)"|>;

check[src_] := Lean`Internal`CheckSource[src, ml, False, None];
body[t_] := "theorem t : True := trivial\n  x " <> t <> " y\n";
col0[t_] := "theorem t : True := trivial\n" <> t <> " x\n";

If[!TrueQ[check[body[""]]["OK"]] || !TrueQ[check[col0["theorem"]]["OK"]],
  Print["FAIL: the probe sources themselves are refused"]; Exit[1]];

If[Environment["MATHILDA_TOKENS_WORDS"] === "1",
  Print[StringRiffle[Sort[DeleteDuplicates[Flatten[StringCases[#,
      RegularExpression["[A-Za-z0-9_]+"]] & /@
    Select[lead, check[body[#]]["Kind"] =!= "lex" &]]]], "\n"]];
  Exit[0]];

missing = {};
Do[
  If[!KeyExistsQ[exempt, t],
    If[TrueQ[check[body[t]]["OK"]], AppendTo[missing, t <> "  (inside a command)"]];
    If[!MemberQ[starters, t] && TrueQ[check[col0[t]]["OK"]],
      AppendTo[missing, t <> "  (at column 0)"]]],
  {t, lead}];
If[nonidx =!= {"nonindexed 0"}, AppendTo[missing, "unindexed command parsers: " <> ToString[nonidx]]];
(* Tokens that end in a word character but are not identifier-shaped: Lean's
   next token starts right after them, possibly inside one of the checker's
   words.  A keyword glued to them must still be refused. *)
odd = Select[toks, StringMatchQ[#, RegularExpression[".*[A-Za-z0-9_]"]] &&
    !StringMatchQ[#, RegularExpression["[A-Za-z_][A-Za-z0-9_'!?]*|[0-9]+"]] &];
Do[If[TrueQ[check[body[t <> w]]["OK"]],
    AppendTo[missing, t <> w <> "  (glued token; add " <> t <> " to LEAN_GLUED_TOKENS)"]],
  {t, odd}, {w, {"theorem", "run_cmd", "1theorem"}}];
Print[If[ml, "Mathlib", "core"], ": ", Length[lead], " command-leading tokens, ",
  Length[toks], " tokens (", Length[odd], " glued), ",
  Length[Keys[exempt]], " exempt (", StringRiffle[Keys[exempt], ", "], ")"];
If[missing === {},
  Print["PASS: every command-leading token is refused by Lean`Internal`CheckSource"]; Exit[0],
  Print["FAIL: not refused (add the word to src/lean_check_tokens.h):"];
  Scan[Print["  ", #] &, missing]; Exit[1]];
