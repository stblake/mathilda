(* ====================================================================== *)
(* LeanBackend.m -- the Lean 4 backend of Certify (loaded by Lean.m).      *)
(*                                                                          *)
(* Discovery ($LeanExecutable, $LeanProject, $LeanVersion), LeanVerify and  *)
(* its source check, user Method scripts, the translation of held           *)
(* statements to Lean text, the per-obligation encoders (leanEncode: the    *)
(* neutral IR of Lean.m's generators -> "Statement", "Proof", "Tier", ...), *)
(* the Tier 1 HasDerivAt walk, chunk / footer / file assembly, runBatch and *)
(* the verdict attribution, the verdict cache, the "Lean" registry entry    *)
(* and the Lean -> Mathilda oracle.  All ::usage strings and messages are   *)
(* in Lean.m.                                                               *)
(*                                                                          *)
(* The C side (src/print_lean.c, src/lean_run.c, src/lean_check.c) supplies *)
(* Lean`Internal`ToLean (the strict LeanForm translator),                   *)
(* Lean`Internal`RunProcess, Lean`Internal`VerifySource (the system()-based *)
(* runner and NDJSON parser) and Lean`Internal`CheckSource.                 *)
(*                                                                          *)
(* The file opens the package itself, with nothing before Begin, so every   *)
(* symbol it creates is in Lean`Private` (no Lean` symbol, nothing in       *)
(* Global`).                                                                *)
(* ====================================================================== *)

BeginPackage["Lean`"];
Begin["`Private`"];

(* ---------------------------------------------------------------------- *)
(* Settings and discovery                                                  *)
(* ---------------------------------------------------------------------- *)

$LeanExecutable = Automatic;

$LeanProject = With[{e = Environment["MATHILDA_LEAN_PROJECT"]},
  If[StringQ[e] && e =!= "", e, None]];

(* The Prelude sits beside this file; $InputFileName is only set while the
   package loads, so the path is captured now (the file is read lazily). *)
$preludeFile = If[StringQ[$InputFileName] && $InputFileName =!= "",
  FileNameJoin[{DirectoryName[$InputFileName], "Prelude.lean"}], $Failed];

(* The kernel-replay checker (Replay.lean, beside this file): every Lean run
   of Certify and LeanVerify compiles the source to an .olean and replays it
   in a fresh lean process (Lean`Internal`VerifySource's fifth argument). *)
$replayFile = If[StringQ[$InputFileName] && $InputFileName =!= "",
  FileNameJoin[{DirectoryName[$InputFileName], "Replay.lean"}], $Failed];

replayScript[] := If[StringQ[$replayFile] && FileExistsQ[$replayFile], $replayFile, None];

(* The core-only prelude (Prelude.lean; no imports, no Mathlib), inlined
   into every generated certificate file and mirrored byte-for-byte as
   lean/Mathilda/Prelude.lean for the Lean-side package (ctest
   lean_prelude_sync).  It carries NO comments or strings: LeanVerify's
   checker skips it as trusted text, and a prelude without any span that Lean
   skips guarantees that Lean's lexer is at the top level where it ends
   (CheckSource refuses a prelude with a string, comment or char literal).
   Its documentation therefore lives here.  Everything is a computable
   definition that the kernel evaluates under `decide +kernel`: recursion is
   structural (an explicit fuel argument where the natural measure is not),
   never well-founded, so the kernel reduces it.
     powModAux fuel b e m acc  square-and-multiply b^e mod m; fuel bounds the
                               loop steps.
     powMod b e m              b^e mod m for every e: the loop halves e and
                               stops at 0, so it runs log2 e + 1 <= e + 1
                               steps and fuel e + 1 never runs out (a fixed
                               fuel silently truncated large exponents; the
                               kernel only peels as many successors off the
                               fuel literal as there are steps).
     prodPow fs                the product of q^k over [(q, k), ...].
     lucasStep known p a fs    one Lucas step: a has order p - 1 mod p, the
                               factor list multiplies out to p - 1, every
                               prime factor q is already known prime.
     prattCheck known steps    a flat Pratt certificate, steps (p, a, factors
                               of p - 1) in dependency order; 2 is known.
     IsPrime n                 the definition: 2 <= n and the only divisors
                               of n are 1 and n (Certify proves IsPrime 2 and
                               Not IsPrime n for n < 2 and composites).
     PrattCertifies p cert     cert is checked and its last step is p.
     Mat, matMul, matVec       dense Rat matrices as row lists; the column
                               count of matMul is taken from B's first row.
     HasShape A r c, HasLength v n
                               exact shapes (matMul / matVec / det use
                               zipWith and getD, which silently truncate or
                               pad ragged input, so every matrix statement
                               states the shapes).
     idMat n, dropCol j A, detAux, det
                               identity; column removal; Laplace expansion
                               along row 0 with the dimension as fuel.
     sqnn, mulnn               Lean.Grind.OrderedRing.sq_nonneg / mul_nonneg. *)
prelude[] := prelude[] = If[StringQ[$preludeFile] && FileExistsQ[$preludeFile],
  ReadString[$preludeFile], ""];

(* The newest ~/.elan toolchain's lean binary, called directly so no
   lean-toolchain file is needed; None when there is none. *)
toolchainVersion[dir_String] := With[{v = StringCases[dir,
     RegularExpression["v(\\d+)\\.(\\d+)\\.(\\d+)"] -> {"$1", "$2", "$3"}]},
  If[v === {}, {0, 0, 0}, ToExpression[First[v]]]];

elanLean[] := Module[{dirs, bins},
  dirs = FileNames["*", FileNameJoin[{$HomeDirectory, ".elan", "toolchains"}]];
  bins = Select[Flatten[{FileNameJoin[{#, "bin", "lean"}],
       FileNameJoin[{#, "bin", "lean.exe"}]} & /@ SortBy[dirs, toolchainVersion]],
    FileExistsQ];
  If[bins === {}, None, Last[bins]]];

(* Automatic discovery, memoized on the first hit. *)
autoLean[] := autoLean[] = Module[{env = Environment["MATHILDA_LEAN"], tc},
  Which[
    StringQ[env] && env =!= "", env,
    (tc = elanLean[]) =!= None, tc,
    True, "lean"]];

leanExe[override_] := Which[
  StringQ[override], override,
  StringQ[$LeanExecutable], $LeanExecutable,
  True, autoLean[]];

leanProject[override_] := Which[
  StringQ[override], override,
  override === None, None,
  StringQ[$LeanProject], $LeanProject,
  True, None];

leanVersionOf[exe_String] := leanVersionOf[exe] = Module[{r},
  r = Quiet[Lean`Internal`RunProcess[{exe, "--version"}, None]];
  If[AssociationQ[r] && r["ExitCode"] === 0 && StringQ[r["StandardOutput"]],
    StringTrim[r["StandardOutput"]], Missing["NotAvailable"]]];

$LeanVersion := leanVersionOf[leanExe[Automatic]];

(* A run that never reached Lean: the shell could not find or execute it. *)
missingLeanQ[r_?AssociationQ] := MemberQ[{126, 127, 9009}, r["ExitCode"]] &&
  AllTrue[r["Messages"], #["Kind"] === "process" &];
missingLeanQ[_] := False;

(* ---------------------------------------------------------------------- *)
(* LeanVerify                                                              *)
(* ---------------------------------------------------------------------- *)

$standardAxioms = {"propext", "Classical.choice", "Quot.sound"};

(* Words a user Method script may not use (Certify; see methodProblem).
   LeanVerify no longer relies on a word list: see checkSource below. *)
$unsafeWords = {"axiom", "native_decide", "implemented_by", "extern", "unsafe",
  "ofReduceBool", "trustCompiler", "debug", "run_cmd", "run_tac", "run_elab", "run_meta",
  "by_elab", "elab", "elab_rules", "macro", "macro_rules", "syntax", "declare_syntax_cat",
  "initialize", "builtin_initialize", "addDeclWithoutChecking", "addDecl", "addAndCompile",
  "skipKernelTC", "command_elab", "term_elab", "tactic_elab", "setEnv", "modifyEnv",
  "withOptions", "native"};

(* The words of a source as LeanVerify's checker splits them: maximal runs of
   [A-Za-z0-9_] ("Lean.Elab.x'" gives Lean, Elab, x), finer than Lean's own
   identifiers. *)
sourceWords[src_String] := DeleteDuplicates[StringCases[src,
   RegularExpression["[A-Za-z0-9_]+"]]];

(* The structural source check (src/lean_check.c).  The guarantee of
   LeanVerify (and of Certify) holds only for sources that cannot run
   user-controlled code while Lean elaborates them: code run at elaboration
   time can forge everything the elaborating process produces, including the
   .olean that the kernel replay reads.  So a source is accepted only if
   (A) it stays inside a minimal lexical surface -- no strings, comments,
   character literals, #, $, %, backquotes, backslashes, @[, French quotes,
   control characters or Unicode outside a short list of mathematical symbols
   and Greek letters -- so that no span exists whose content Lean would skip
   and the checker need not reproduce Lean's lexer; (B) none of its
   [A-Za-z0-9_] words is a command-leading word of Lean's parser (outside the
   allowed command starters at column 0; the list is derived from Lean's
   parser tables, src/lean_check_tokens.h) or an escape word (by_elab,
   native_decide, sorry, eval, IO, ...); and (C) it fits an ALLOWLIST grammar
   -- no imports (core Lean only), or only `import Mathlib[.Path]` with a
   project; only theorem / lemma (Mathlib) / example (no open / namespace:
   they activate scoped syntax) / section / end / variable / universe / numeric or Boolean set_option
   commands at column 0.  The trusted prelude (Prelude.lean, byte-identical,
   itself free of strings and comments) is skipped.  With certify -> True the
   generator auxiliaries `def mathilda_cert_<k>_<x>` and the footer lines
   `#print axioms mathilda_cert_<k>_footer` are accepted too.  Gives the C
   checker's association. *)
checkSource[src_String, mathlib_, certify_] :=
  Lean`Internal`CheckSource[src, TrueQ[mathlib], TrueQ[certify],
    With[{p = prelude[]}, If[StringQ[p] && p =!= "", p, None]]];

(* The message of each refusal kind. *)
refusalMessage["import", w_, _] := Message[LeanVerify::import, w];
refusalMessage["cmd", w_, _] := Message[LeanVerify::cmd, w];
refusalMessage["token", w_, _] := Message[LeanVerify::unsafe, w];
refusalMessage[_, w_, l_] := Message[LeanVerify::lex, w, l];

(* src (examples already renamed by the checker) followed by a trusted footer:
   one `#print axioms` line per declared name.  Returns {source, <|name ->
   footer line|>}; only a report printed on its own footer line counts (a
   report elsewhere could be forged). *)
withAxiomReports[src_String, names_List] := Module[{base, l0},
  base = If[StringEndsQ[src, "\n"] || src === "", src, src <> "\n"];
  l0 = StringCount[base, "\n"] + 1;
  {StringJoin[base, ("#print axioms " <> # <> "\n") & /@ names],
   AssociationThread[names -> Range[l0, l0 + Length[names] - 1]]}];

anchoredAxioms[r_, name_, line_] := With[{rep = Select[Lookup[r, "AxiomReports", {}],
     #["Theorem"] === name && #["Line"] === line &]},
  If[Length[rep] === 1, rep[[1]]["Axioms"], Missing["NotReported"]]];

(* Option names are strings (v0.6523): OptionValue matches a name by its
   text, so the user's Engine -> v (a Global` symbol) and "Engine" -> v both
   work, and the package creates no option-name symbols of its own -- a
   symbol in Lean` would sit on $ContextPath and capture the user's (or a
   later package's) Engine, Tier, Domain, ... .  Method, TimeConstraint and
   Assumptions are System` symbols. *)
Options[LeanVerify] = {"LeanExecutable" -> Automatic, "LeanProject" -> Automatic,
  "KeepFiles" -> False};

LeanVerify[src_String, opts:OptionsPattern[]] := Module[{exe, proj, full, lines, r, chk,
    ok},
  If[optionNameProblem[LeanVerify, {opts}] =!= None, Return[$Failed, Module]];
  exe = leanExe[OptionValue["LeanExecutable"]];
  proj = leanProject[OptionValue["LeanProject"]];
  chk = checkSource[src, StringQ[proj], False];
  If[!AssociationQ[chk], Return[$Failed, Module]];
  If[!TrueQ[chk["OK"]],
    refusalMessage[chk["Kind"], chk["What"], chk["Line"]];
    Return[<|"Success" -> False, "ExitCode" -> Missing["NotRun"],
      "Messages" -> {<|"Severity" -> "error", "Kind" -> "mathilda." <> chk["Kind"],
        "Line" -> chk["Line"], "Text" -> "not run: the source uses " <> chk["What"]|>},
      "Axioms" -> <||>, "AxiomReports" -> {}, "Timing" -> 0., "Command" -> "",
      "Rejected" -> {chk["What"]}|>, Module]];
  {full, lines} = withAxiomReports[chk["Source"], chk["Theorems"]];
  r = Lean`Internal`VerifySource[full, exe, proj, TrueQ[OptionValue["KeepFiles"]],
    replayScript[]];
  If[!AssociationQ[r], Return[$Failed, Module]];
  If[missingLeanQ[r], Message[LeanVerify::nolean]];
  (* Success: Lean's own run was clean, the kernel replay of the compiled
     module accepted every declaration with standard axioms, and every
     theorem / example the SOURCE TEXT declares (names from the checker, not
     from Lean's output) is among the replayed declarations with standard
     axioms and has a standard report on its own footer line -- so a module
     that drops or renames a declaration fails. *)
  ok = TrueQ[r["Success"]] && r["Replay"] === True && AllTrue[Keys[lines], Function[nm,
    With[{ax = anchoredAxioms[r, nm, lines[nm]], rax = Lookup[r["ReplayAxioms"], nm]},
      ListQ[ax] && SubsetQ[$standardAxioms, ax] &&
        ListQ[rax] && SubsetQ[$standardAxioms, rax]]]];
  Join[<|"Success" -> ok|>, KeyDrop[r, "Success"]]];

(* ---------------------------------------------------------------------- *)
(* Translation and source assembly                                         *)
(* ---------------------------------------------------------------------- *)


(* Strict translation of a held statement; $Failed (quietly) when refused.
   closed -> True (every theorem statement) refuses free symbols. *)
toLean[held_HoldForm, ty_, coreOnly_, pass_List:{}, closed_:False] :=
  Quiet[Lean`Internal`ToLean[held, ty, coreOnly, pass, TrueQ[closed]]];

heartbeats[t_?NumericQ] /; t > 0 := Round[40000 t];
heartbeats[_] := 0;

(* ---- User-supplied Method strings ----------------------------------------
   A Method string is pasted into the generated file as a tactic block, so it
   must not be able to leave that block or change how later text is read.
   It must satisfy LeanVerify's own rules (A) and (B) above: methodProblem
   runs the C checker on `theorem t : True := by` followed by the script,
   indented by two spaces, so the script inherits the minimal lexical surface
   (no strings, comments, #, $, %, @[, backquotes, ...) and the refusal of
   every command-leading and escape word.  On top of that the words below
   (Lean, IO, axiom, ...) are refused, and every line is indented by two
   spaces, so no line reaches column 0; Certify runs every item with a user
   Method in a Lean file of its own (certifyItems), after which the trusted
   footer re-checks the theorem's statement and axioms. *)
$methodBannedWords = Join[$unsafeWords, {"theorem", "lemma", "def", "abbrev", "instance",
  "example", "structure", "class", "inductive", "notation", "infix", "infixl", "infixr",
  "prefix", "postfix", "set_option", "attribute", "open", "export", "namespace",
  "section", "end", "universe", "variable", "import", "mutual", "deriving", "opaque",
  "partial", "private", "protected", "noncomputable", "omit", "include", "scoped",
  "local", "sorry", "admit", "Lean", "IO"}];
$methodBannedText = {"\[LeftGuillemet]", "\[RightGuillemet]", "`", "$", "#", "@[", "--",
  "/-", "-/", "\"", "\t", "\r", "%", "\\"};

methodProblem[m_String] := Module[{w, t, chk},
  If[StringTrim[m] === "", Return["the script is empty", Module]];
  t = Select[$methodBannedText, StringContainsQ[m, #] &];
  If[t =!= {}, Return["it contains " <> ToString[First[t], InputForm], Module]];
  If[StringContainsQ[m, RegularExpression["[\\x00-\\x09\\x0B-\\x1F\\x7F]"]],
    Return["it contains a control character", Module]];
  w = Intersection[sourceWords[m], $methodBannedWords];
  If[w =!= {}, Return["it uses the keyword " <> First[w], Module]];
  chk = Lean`Internal`CheckSource["theorem mathilda_method : True := by\n" <>
    StringRiffle[("  " <> #) & /@ StringSplit[m, "\n"], "\n"] <> "\n", False, False, None];
  If[!TrueQ[chk["OK"]], Return["it uses " <> ToString[chk["What"]] <>
    " (LeanVerify's source rules apply to Method scripts)", Module]];
  None];
methodProblem[m_] := If[MemberQ[{Automatic, "grind", "decide", "omega"}, m], None,
  "Method must be Automatic, \"grind\", \"decide\", \"omega\" or a tactic script string"];

(* A user script: each line indented by two spaces under `by`. *)
userMethodQ[m_] := StringQ[m] && !MemberQ[{"grind", "decide", "omega"}, m];

proofScript[method_] := Which[
  method === Automatic || method === "grind", "by\n  intros\n  grind",
  method === "decide", "by\n  decide",
  method === "omega", "by\n  intros\n  omega",
  StringQ[method], "by\n" <> StringRiffle[("  " <> #) & /@ StringSplit[method, "\n"], "\n"],
  True, "by\n  intros\n  grind"];

(* One theorem chunk: <|"Name", "Text", "Statement"|>.  A generator may also
   supply a "Preamble" (auxiliary defs / lemmas placed before the theorem).
   The statement, proof and preamble name their auxiliary declarations
   @NAME@_..., replaced here by the theorem's name.  The chunk's
   `#print axioms` is not in the chunk: it is in the trusted footer
   (footerText), after all user-influenced text. *)
theoremChunk[name_String, stmt_String, proof_String, hb_Integer, pre_String:""] :=
  With[{sp = splitGeneric[stmt]}, <|
  "Name" -> name,
  "Statement" -> StringReplace[stmt, "@NAME@" -> name],
  "Text" -> StringReplace[StringJoin[pre,
    If[hb > 0, "set_option maxHeartbeats " <> ToString[hb] <> " in\n", ""],
    "theorem ", name, sp[[1]], " :\n    ", sp[[2]], " := ", proof, "\n"],
    "@NAME@" -> name]|>];

(* The footer lines of one chunk: the statement is re-elaborated and checked
   against the theorem (so a chunk cannot pass by declaring something else
   under that name), then its axioms are reported. *)
footerName[c_] := c["Name"] <> "_footer";
footerText[c_] := StringJoin["theorem ", footerName[c], " :\n    ", c["Statement"], " :=\n  ",
  c["Name"], "\n", "#print axioms ", footerName[c], "\n"];

fileHeader[tier_] := StringJoin[
  If[tier === 1, "import Mathlib\nset_option linter.unusedTactic false\n\
set_option linter.unreachableTactic false\n", ""],
  "set_option autoImplicit false\n\n"];

(* The file for a list of chunks plus, per chunk, its line spans {{first,
   last} of the chunk, {first, last} of its footer} and the line of its
   footer's `#print axioms`. *)
assembleFile[chunks_List, tier_] := Module[{head, line, spans = {}, body = {}, foot = {},
    fspans = {}, axl = {}},
  head = fileHeader[tier] <> prelude[] <> "\n";
  line = StringCount[head, "\n"] + 1;
  Do[
    With[{t = c["Text"]},
      AppendTo[spans, {line, line + StringCount[t, "\n"] - 1}];
      AppendTo[body, t <> "\n"];
      line += StringCount[t, "\n"] + 1],
    {c, chunks}];
  Do[
    With[{t = footerText[c]},
      AppendTo[fspans, {line, line + StringCount[t, "\n"] - 1}];
      AppendTo[axl, line + StringCount[t, "\n"] - 1];
      AppendTo[foot, t];
      line += StringCount[t, "\n"]],
    {c, chunks}];
  {StringJoin[head, body, foot], Transpose[{Transpose[{spans, fspans}], axl}]}];

(* ---------------------------------------------------------------------- *)
(* Running Lean and attributing verdicts (the Verified contract)           *)
(* ---------------------------------------------------------------------- *)

badMessageQ[m_] := m["Severity"] === "error" || m["Kind"] === "hasSorry";
inSpanQ[m_, {a_, b_}] := a <= m["Line"] <= b;
inSpansQ[m_, sp_List] := AnyTrue[sp, inSpanQ[m, #] &];

(* Verdict of the theorem `name`, given its layout {spans, axiom line} and
   the spans of every chunk in the file.  True needs all of:
     - no error / sorry inside its own spans (its chunk and its footer);
     - no error / sorry outside every chunk's spans (prelude, header);
     - the `#print axioms` report printed on its own footer line (a report
       anywhere else is ignored), naming only the standard axioms;
     - Lean exited with code 0, or every error it reported lies in some
       chunk (then the non-zero exit is explained by other chunks);
     - the kernel replay of the compiled file (Replay.lean, a fresh lean
       process that runs none of the file's code) accepted every
       declaration, and found only standard axioms under the footer theorem.
   A file with an error in some chunk is not compiled, so it cannot be
   replayed: runBatch then re-runs the clean chunks on their own. *)
attribute[r_, name_, {span_, axLine_}, spans_] := Module[{bad, own, stray, ax, rax, logs},
  bad = Select[r["Messages"], badMessageQ];
  own = Select[bad, inSpansQ[#, span] &];
  stray = Select[bad, Function[m, NoneTrue[spans, inSpansQ[m, #] &]]];
  ax = anchoredAxioms[r, name <> "_footer", axLine];
  rax = If[r["Replay"] === True, Lookup[r["ReplayAxioms"], name <> "_footer"], Missing[]];
  logs = Function[m, <|"Severity" -> m["Severity"], "Line" -> m["Line"],
      "Text" -> m["Text"]|>] /@
    Select[r["Messages"], Function[m, inSpansQ[m, span] ||
        (badMessageQ[m] && NoneTrue[spans, inSpansQ[m, #] &])]];
  <|"Verified" -> (own === {} && stray === {} && ListQ[ax] &&
      SubsetQ[$standardAxioms, ax] && (r["ExitCode"] === 0 || bad =!= {}) &&
      ListQ[rax] && SubsetQ[$standardAxioms, rax]),
    "Axioms" -> If[ListQ[ax], ax, {}], "Logs" -> logs|>];

(* Run one file of chunks; a list of per-chunk verdict associations.  When
   Lean rejected some chunks, the file was not compiled and so not replayed;
   the chunks with no error of their own are then run again in a file of
   their own (their verdict needs the replay), the others stay False. *)
runBatch[chunks_List, tier_, exe_, proj_] := Module[{src, lay, r, res, clean, chk},
  {src, lay} = assembleFile[chunks, tier];
  (* Internal assertion: the generated file passes the same allowlist check as
     a LeanVerify source (with the trusted prelude, the Tier 1 Mathlib import
     and the generator auxiliaries def mathilda_cert_<k>_<x>).  The only
     user-influenced text in it is a Method script (deny-listed) and the
     strictly translated statements. *)
  chk = checkSource[src, tier === 1, True];
  If[!TrueQ[chk["OK"]],
    Message[Certify::unsafe, chk["What"]];
    Return[Table[<|"Verified" -> False, "Axioms" -> {},
      "Logs" -> {<|"Severity" -> "error", "Line" -> chk["Line"],
        "Text" -> "not run: the source check refused " <> ToString[chk["What"]]|>},
      "Timing" -> 0., "Source" -> src|>, {Length[chunks]}], Module]];
  r = Lean`Internal`VerifySource[src, exe, If[tier === 1, proj, None], False, replayScript[]];
  If[AssociationQ[r] && !missingLeanQ[r] && r["Replay"] =!= True && r["ExitCode"] =!= 0 &&
      Length[chunks] > 1,
    res = MapThread[attribute[r, #1["Name"], #2, lay[[All, 1]]] &, {chunks, lay}];
    clean = Select[Range[Length[chunks]], Function[k, NoneTrue[res[[k]]["Logs"],
      badMessageQ[Append[#, "Kind" -> ""]] &]]];
    If[0 < Length[clean] < Length[chunks],
      Return[Module[{again = AssociationThread[clean -> runBatch[chunks[[clean]], tier, exe, proj]]},
        Table[If[KeyExistsQ[again, k], again[k],
          Join[res[[k]], <|"Verified" -> False, "Timing" -> r["Timing"], "Source" -> src|>]],
          {k, Length[chunks]}]], Module]]];
  Which[
    !AssociationQ[r],
      Table[<|"Verified" -> Missing["NotAvailable"], "Axioms" -> {}, "Logs" -> {},
        "Timing" -> 0., "Source" -> src|>, {Length[chunks]}],
    missingLeanQ[r],
      Message[LeanVerify::nolean];
      Table[<|"Verified" -> Missing["NotAvailable"], "Axioms" -> {},
        "Logs" -> (KeyDrop[#, "Kind"] & /@ r["Messages"]), "Timing" -> r["Timing"],
        "Source" -> src|>, {Length[chunks]}],
    True,
      MapThread[Join[attribute[r, #1["Name"], #2, lay[[All, 1]]],
          <|"Timing" -> r["Timing"] + Lookup[r, "ReplayTiming", 0.], "Source" -> src|>] &, {chunks, lay}]]];

(* The generic form (Lean.m, identityGen): a statement translated over Rat
   is generalised, the word Rat becoming the type variable K bound by
   $charZeroBinders / $orderedBinders (core Lean's grind classes).  The
   theorem takes them as parameters; the footer states the full
   forall {K : Type} [...], ... form. *)
$charZeroBinders = "{K : Type} [Lean.Grind.Field K] [Lean.Grind.IsCharP K 0]";
$orderedBinders = "{K : Type} [Lean.Grind.Field K] [LE K] [LT K] [Std.LawfulOrderLT K] \
[Std.IsLinearOrder K] [Lean.Grind.OrderedRing K]";
genericBinders["CharZero"] = $charZeroBinders; genericBinders["Ordered"] = $orderedBinders;

(* Rat (as a whole word, not Rat.f) -> K, in a statement or a proof script. *)
genericText[s_String] := StringReplace[s,
  RegularExpression["(?<![A-Za-z0-9_.'])Rat(?![A-Za-z0-9_.'])"] -> "K"];
genericText[other_] := other;
genericStmt[s_String, kind_String] := "\[ForAll] " <> genericBinders[kind] <> ", " <> genericText[s];
genericStmt[other_, _] := other;

(* {binders, body} of a generic statement; {"", s} for any other. *)
splitGeneric[s_String] := Module[{b = SelectFirst[{$charZeroBinders, $orderedBinders},
     StringStartsQ[s, "\[ForAll] " <> # <> ", "] &, None]},
  If[b === None, {"", s}, {" " <> b, StringDrop[s, StringLength["\[ForAll] " <> b <> ", "]]}]];

(* Prelude functions in statements are written with private placeholder
   heads and renamed after translation.  A symbol in a context on
   $ContextPath (Mathilda`matMul) would capture the user's short name. *)
$preludeHeads = {"Lean.Private.matMul" -> "Mathilda.matMul",
  "Lean.Private.matVec" -> "Mathilda.matVec", "Lean.Private.idMat" -> "Mathilda.idMat",
  "Lean.Private.det" -> "Mathilda.det", "Lean.Private.hasShape" -> "Mathilda.HasShape",
  "Lean.Private.hasLength" -> "Mathilda.HasLength"};

(* Private generator variables (the surd symbols s1$n, s2$n, ...) print as
   <<Lean`Private`s1$n>>; they are renamed to their bare names s1, s2, ... *)
leanRename[s_String] := StringReplace[StringReplace[StringReplace[s, $preludeHeads],
  {"\[LeftGuillemet]Lean`Private`" ~~ n:(LetterCharacter ~~ DigitCharacter..) ~~
     "$" ~~ DigitCharacter.. ~~ "\[RightGuillemet]" :> n,
   "\[LeftGuillemet]Lean`Private`" ~~ n:(LetterCharacter ~~ DigitCharacter..) ~~
     "\[RightGuillemet]" :> n}],
  (* Lean`Private`ratCast[m1]: the Nat binder m1 used as a rational value. *)
  {"Lean.Private.ratCast " ~~ n:(LetterCharacter ~~ DigitCharacter..) :>
    "(" <> n <> " : Rat)", "Lean.Private.fib" -> "@NAME@_F"}];
leanRename[other_] := other;

(* Strict Tier 0 translation of a held (sub)term over type ty; its free
   symbols are the caller's business (proof-script pieces under intro). *)
leanStmt[held_HoldForm, ty_String, extra_List:{}] :=
  leanRename[toLean[held, ty, True, Join[Keys[$preludeHeads], extra]]];

(* The same for a whole theorem statement: every symbol must be bound by a
   quantifier of the statement (a free `default` would be Lean's
   Inhabited.default, a free parameter an unknown identifier). *)
leanThm[held_HoldForm, ty_String, extra_List:{}] :=
  leanRename[toLean[held, ty, True, Join[Keys[$preludeHeads], extra], True]];

(* Lean name of a binder variable. *)
leanVar[v_Symbol] := With[{s = leanStmt[HoldForm[v], "Rat"]}, If[StringQ[s], s, "_"]];

introLine[{}] := "";
introLine[vars_List] := "\n  intro " <> StringRiffle[leanVar /@ vars, " "];

ratText[e_] := leanStmt[HoldForm[e], "Rat"];

(* Lean `have` lines for the square hints; gnames maps each multiplier g to
   the name of its hypothesis 0 <= g. *)
sqHints[terms_, gnames_] := Table[With[{t = terms[[k]]},
    "\n  have mth_s" <> ToString[k] <> " : 0 \[LessEqual] (" <> ratText[t[[2]]] <> ") ^ 2" <>
      If[t[[3]] === 1, " := Mathilda.sqnn",
        " * (" <> ratText[t[[3]]] <> ") := Mathilda.mulnn Mathilda.sqnn " <> (t[[3]] /. gnames)]],
  {k, Length[terms]}];

sqSum[terms_] := StringRiffle[("(" <> ratText[#[[1]]] <> ") * (" <> ratText[#[[2]]] <> ") ^ 2" <>
    If[#[[3]] === 1, "", " * (" <> ratText[#[[3]]] <> ")"]) & /@ terms, " + "];

(* The remaining conjuncts (feasibility, value) are closed rational facts. *)
nestAnd[1] := "(by grind)";
nestAnd[n_Integer] := "(And.intro (by grind) " <> nestAnd[n - 1] <> ")";

prattPair[{q_, k_}] := "(" <> ToString[q] <> ", " <> ToString[k] <> ")";
prattStep[{p_, a_, fs_}] := "(" <> ToString[p] <> ", " <> ToString[a] <> ", [" <>
  StringRiffle[prattPair /@ fs, ", "] <> "])";
prattText[cert_List] := "[" <> StringRiffle[prattStep /@ cert, ", "] <> "]";

$lucasStep = "Lucas's primality theorem (Tier 0 checks the Pratt certificate: Lean evaluates \
Mathilda.prattCheck, it does not prove Lucas's theorem; 2 is the base of the scheme)";

nestAndTerms[{t_}] := t;
nestAndTerms[{t_, r__}] := "And.intro (" <> t <> ") (" <> nestAndTerms[{r}] <> ")";

factorPow[{p_, 1}] := ToString[p];
factorPow[{p_, k_}] := ToString[p] <> " ^ " <> ToString[k];

(* c^(alpha v + beta) (c rational, alpha, beta integers) with v -> B + sh as
   (c^alpha)^B * c^(alpha sh + beta): the exponent is the Nat binder B only. *)
natPower[c_, ex_, v_, b_, sh_] := Module[{al = Coefficient[Expand[ex], v, 1], be},
  be = Expand[ex - al v];
  If[c =!= 0 && IntegerQ[al] && IntegerQ[be] && FreeQ[be, v],
    iProd[iPower[c^al, b], c^(al sh + be)], $Failed]];

(* e with the summation variable v replaced by the Nat binder b shifted by
   sh: exponents become b (natPower), every other v the rational (b + sh). *)
natForm[e_, v_, b_, sh_] := Module[{r},
  r = e /. Power[c_?ratNumberQ, ex_] /; !FreeQ[ex, v] :> natPower[c, ex, v, b, sh];
  r = r /. Fibonacci[u_] :> natFib[Expand[u - v] + sh, b];
  If[!FreeQ[r, $Failed], Return[$Failed, Module]];
  r /. v -> If[sh === 0, ratCast[b], iPlus[ratCast[b], sh]]];

(* Fibonacci[v + j] as the source's _F applied to the Nat binder plus the
   natural shift j + sh (placeholder fib, renamed by leanRename). *)
natFib[j_Integer, b_] /; j >= 0 := If[j === 0, fib[b], fib[iPlus[b, j]]];
natFib[_, _] := $Failed;

$sumBinder = Lean`Private`m1;
$termBinder = Lean`Private`k1;
$sumPass = {"Lean.Private.ratCast", "Lean.Private.m1", "Lean.Private.k1", "Lean.Private.fib"};

(* ---------------------------------------------------------------------- *)
(* Row 17: Tier 1 antiderivative (Mathlib) -- HasDerivAt F f x, the proof   *)
(* term built structurally from F's tree (plan section 3.3) and closed by   *)
(* congr_deriv + field_simp / ring.  The function body is written by the    *)
(* same walk, so it matches the proof term's shape exactly.                 *)
(* ---------------------------------------------------------------------- *)

(* Lean text of a constant (free of x) as a Real; $Failed if untranslatable. *)
realConst[c_] := With[{s = toLean[HoldForm[c], "Real", False, {}, True]},
  If[StringQ[s], "(" <> s <> " : Real)", $Failed]];

(* {body, proof term} for e (in x); hypotheses d != 0 are appended to the
   private list $hdHyps as {name, body of d, d}.  $Failed outside the table. *)
hd[e_, x_] /; FreeQ[e, x] := With[{c = realConst[e]},
  If[c === $Failed, $Failed, {c, "(hasDerivAt_const " <> $hdVar <> " " <> c <> ")"}]];
hd[x_Symbol, x_Symbol] := {$hdVar, "(hasDerivAt_id' " <> $hdVar <> ")"};
hd[e_Plus, x_] := Module[{ps = hd[#, x] & /@ (List @@ e)},
  If[MemberQ[ps, $Failed], $Failed,
    Fold[{"(" <> #1[[1]] <> " + " <> #2[[1]] <> ")", "(" <> #1[[2]] <> ".add " <> #2[[2]] <> ")"} &,
      First[ps], Rest[ps]]]];
hd[e_Times, x_] := Module[{fs = List @@ e, cs, num, den, n, d, pn, pd, h},
  cs = Select[fs, FreeQ[#, x] &];
  den = Cases[fs, Power[b_, k_Integer?Negative] /; !FreeQ[b, x] :> b^(-k)];
  num = Select[fs, !FreeQ[#, x] && !MatchQ[#, Power[b_, _Integer?Negative] /; !FreeQ[b, x]] &];
  Which[
    den =!= {},
      pn = hd[Times @@ Join[cs, num], x]; pd = hd[Times @@ den, x];
      If[pn === $Failed || pd === $Failed, Return[$Failed, Module]];
      h = hdHyp[pd[[1]], Times @@ den];
      {"(" <> pn[[1]] <> " / " <> pd[[1]] <> ")", "(" <> pn[[2]] <> ".div " <> pd[[2]] <> " " <> h <> ")"},
    cs =!= {},
      n = realConst[Times @@ cs]; pn = hd[Times @@ num, x];
      If[n === $Failed || pn === $Failed, Return[$Failed, Module]];
      {"(" <> n <> " * " <> pn[[1]] <> ")", "(" <> pn[[2]] <> ".const_mul " <> n <> ")"},
    True,
      pn = hd[#, x] & /@ num;
      If[MemberQ[pn, $Failed], Return[$Failed, Module]];
      Fold[{"(" <> #1[[1]] <> " * " <> #2[[1]] <> ")", "(" <> #1[[2]] <> ".mul " <> #2[[2]] <> ")"} &,
        First[pn], Rest[pn]]]];
hd[Power[x_Symbol, k_Integer?Positive], x_Symbol] :=
  {"(" <> $hdVar <> " ^ " <> ToString[k] <> ")",
   "(hasDerivAt_pow " <> ToString[k] <> " " <> $hdVar <> ")"};
hd[Power[u_, k_Integer?Positive], x_] := With[{p = hd[u, x]},
  If[p === $Failed, $Failed,
    {"(" <> p[[1]] <> " ^ " <> ToString[k] <> ")", "(" <> p[[2]] <> ".pow " <> ToString[k] <> ")"}]];
hd[Power[u_, k_Integer?Negative], x_] := Module[{p = hd[u^(-k), x], h},
  If[p === $Failed, Return[$Failed, Module]];
  h = hdHyp[p[[1]], u^(-k)];
  {"(1 / " <> p[[1]] <> ")", "((hasDerivAt_const " <> $hdVar <> " (1 : Real)).div " <> p[[2]] <> " " <> h <> ")"}];
hd[Log[u_], x_] := Module[{p = hd[u, x], h},
  If[p === $Failed, Return[$Failed, Module]];
  h = hdHyp[p[[1]], u];
  {"Real.log " <> p[[1]], "(" <> p[[2]] <> ".log " <> h <> ")"}];
hd[ArcTan[u_], x_] := With[{p = hd[u, x]},
  If[p === $Failed, $Failed, {"Real.arctan " <> p[[1]], "(" <> p[[2]] <> ".arctan)"}]];
(* Mathlib has no Real.artanh: artanh u = (log (1 + u) - log (1 - u)) / 2. *)
hd[ArcTanh[u_], x_] := Module[{p = hd[u, x], h1, h2, a, b},
  If[p === $Failed, Return[$Failed, Module]];
  a = "(1 + " <> p[[1]] <> ")"; b = "(1 - " <> p[[1]] <> ")";
  h1 = hdHyp[a, 1 + u]; h2 = hdHyp[b, 1 - u];
  {"((Real.log " <> a <> " - Real.log " <> b <> ") / 2)",
   "((((" <> p[[2]] <> ".const_add 1).log " <> h1 <> ").sub ((" <> p[[2]] <>
     ".const_sub 1).log " <> h2 <> ")).div_const 2)"}];
hd[_, _] := $Failed;

$hdHyps = {};
$hdVar = "x";

(* The hypothesis d != 0 for a polynomial d, by its Lean text; its name. *)
realHyp[d_] := With[{s = toLean[With[{dd = d}, HoldForm[dd]], "Real", False, {}]},
  If[StringQ[s], hdHyp["(" <> s <> ")", d], $Failed]];
hdHyp[body_String, d_] := Module[{old = SelectFirst[$hdHyps,
    #[[2]] === body || Expand[#[[3]] - d] === 0 &, None]},
  If[old =!= None, old[[1]],
    With[{nm = "mth_h" <> ToString[Length[$hdHyps] + 1]},
      AppendTo[$hdHyps, {nm, body, d}]; nm]]];

(* Real.log semantics for the witness check: Real.log u = Log[Abs[u]],
   Real.log 0 = 0. *)
realLogAt[u_?NumericQ] := If[PossibleZeroQ[u], 0, Log[Abs[u]]];
$witnessPoints = {3, 7/2, 11, 1/3, 5/7, -3, -7/2, -11, -1/3, 17/5, 23, -13/4};
(* A test point x0 at which every d (in x) is nonzero under Lean's Real
   semantics, or $Failed.  The value is decided numerically at 40 digits
   (|d(x0)| > 10^-20); for a hypothesis outside the polynomials in x this is
   recorded as a trusted step (a nonzero polynomial is never identically 0). *)
hypWitness[ds_List, x_Symbol] := SelectFirst[$witnessPoints, Function[x0,
    AllTrue[ds, With[{v = Quiet[N[(# /. Log -> realLogAt) /. x -> x0, 40]]},
      NumericQ[v] && Im[v] == 0 && Abs[v] > 10^-20] &]], $Failed];
$hdWitness = None;

tier1AntiderivGen[input_, F_, f_, x_Symbol, o_] := Module[{p, ftxt, b, hb, hyps, stmt, held,
    xs, proof},
  $hdHyps = {}; $hdWitness = None;
  $hdVar = xs = leanVar[x];
  p = hd[F, x];
  ftxt = toLean[HoldForm[f], "Real", False, {}];
  If[p === $Failed || !StringQ[ftxt] || !FreeQ[F, Complex | Root | RootSum],
    Return[$Failed, Module]];
  (* field_simp also needs the integrand's denominator to be nonzero. *)
  b = Denominator[Together[f]];
  hb = If[FreeQ[b, x], None, realHyp[b]];
  If[hb === $Failed, Return[$Failed, Module]];
  (* Each irreducible factor of a hypothesis is nonzero as well (implied, so
     the statement is not weakened); field_simp needs them one by one. *)
  Do[Do[If[!FreeQ[fc, x] && realHyp[fc] === $Failed, Return[$Failed, Module]],
      {fc, FactorList[d][[All, 1]]}],
    {d, $hdHyps[[All, 3]]}];
  hyps = $hdHyps;
  (* The statement is vacuous if some hypothesis d != 0 fails at every real
     x: Log[x^2] - 2 Log[x] is never 0 over C for x < 0 but is identically 0
     under Real.log (Real.log x = log |x|), so "∀ x, d ≠ 0 → ..." would hold
     for any F.  Require a witness point where every hypothesis is nonzero
     under Lean's semantics. *)
  If[hyps =!= {},
    With[{w = hypWitness[hyps[[All, 3]], x]},
      If[w === $Failed,
        Message[Certify::vacuous, input,
          SelectFirst[hyps[[All, 3]], hypWitness[{#}, x] === $Failed &, hyps[[1, 3]]]];
        Return["vacuous", Module]];
      $hdWitness = w]];
  stmt = StringJoin["\[ForAll] ", xs, " : Real, ",
    (#[[2]] <> " \[NotEqual] 0 \[RightArrow] ") & /@ hyps,
    "HasDerivAt (fun ", xs, " : Real => ", p[[1]], ") ((", ftxt, " : Real)) ", xs];
  proof = StringJoin["by\n  intro ", xs, (" " <> #[[1]]) & /@ hyps,
    "\n  have mth_d := ", p[[2]],
    "\n  refine mth_d.congr_deriv ?_",
    (* HasDerivAt.add etc. state the function as f + g: apply it pointwise. *)
    "\n  try simp only [Pi.add_apply, Pi.sub_apply, Pi.mul_apply, Pi.div_apply]",
    If[hb === None, "", "\n  try rw [eq_div_iff " <> hb <> "]"],
    "\n  try field_simp",
    (* ring1 (ring would fall back to ring_nf without failing), then ring_nf +
       Real.sq_sqrt for surd coefficients such as 1 / Real.sqrt 2. *)
    "\n  all_goals first\n  | ring1\n  | (ring_nf; norm_num [Real.sq_sqrt]; try ring_nf)"];
  held = With[{FF = F, ff = f, xx = x, cond = And @@ (Unequal[#[[3]], 0] & /@ hyps)},
    HoldForm[ForAll[xx, cond, D[FF, xx] == ff]]];
  <|"Theorem" -> held, "Statement" -> stmt, "Proof" -> proof,
    "Certificate" -> <|"Form" -> "HasDerivAt", "Variable" -> x, "ProofTerm" -> p[[2]],
      "HypothesisWitness" -> $hdWitness|>,
    "TrustedSteps" -> If[$hdWitness === None || AllTrue[hyps[[All, 3]], PolynomialQ[#, x] &], {},
      {"the hypotheses are satisfiable: each is nonzero at " <>
      ToString[x] <> " = " <> ToString[$hdWitness, InputForm] <> " under Real.log semantics \
(checked by Mathilda numerically, 40 digits; not by Lean)"}],
    "Hypotheses" -> (With[{d = #[[3]]}, HoldForm[d != 0]] & /@ hyps), "Domain" -> Reals,
    "Claim" -> {D[F, x], f}, "Tier" -> 1|>];

(* Session cache of Lean verdicts, keyed by the theorem text under a fixed
   name (so a theorem is found whatever its position in a batch), the tier,
   the executable and the project.  Only definitive verdicts are stored. *)
$certCache = <||>;
(* The key is the full text itself (not a hash), so two different theorems
   can never share a verdict. *)
cacheKey[g_, hb_, exe_, proj_] := {
  theoremChunk["mathilda_cert", g["Statement"], g["Proof"], hb,
    Lookup[g, "Preamble", ""]]["Text"],
  g["TierUsed"], exe, If[g["TierUsed"] === 1, proj, None]};

(* ---------------------------------------------------------------------- *)
(* Lean -> Mathilda oracle (server side of the lean/ package tactics)      *)
(* ---------------------------------------------------------------------- *)

(* Every answer is plain data as a String (or "MATHILDA-FAIL: reason"):
   numerals, the variables v1, v2, ... that the Lean serializer sent,
   + - * / ^ (natural-number exponents), lists [..] and tuples (..).  The
   Lean side (lean/Mathilda/Answer.lean) parses it with its own small parser
   and never elaborates it as Lean syntax; anything else is refused.
   Coefficients carry no type: Lean builds them at the goal's type.  Nothing
   here is trusted; the tactic re-checks the certificate in the Lean kernel. *)

(* The entry point is Lean`Internal`Oracle, called by its full name from the
   lean/ tactics (Mathilda.Tactics); Lean`Internal` is not on $ContextPath, so
   it does not capture a user's (or a later package's) symbol Oracle. *)
Lean`Internal`Oracle::usage =
  "Lean`Internal`Oracle[kind, args...] is the server side of the Lean 4 tactics in the lean/ \
package (mathilda_linear_combination, mathilda_sos, mathilda_prime, #mathilda). It answers \
with plain data as a String (numerals, v1, v2, ..., + - * / ^, lists and tuples), which the \
Lean side parses itself and never elaborates as Lean syntax: kind \"LinearCombination\" \
({a, b}, {{l1, r1}, ...}, vars) gives cofactors c with a - b == Sum c_i (l_i - r_i); \"SOS\" \
({p, q}, {g1, ...}, vars) gives (d, l, i) triples with q - p == Sum d l^2 g_i (g_0 = 1); \
\"Pratt\" (p) gives a flat Pratt certificate; \"Eval\" (\"input\") gives the LeanForm of the \
evaluated input. A failure gives \"MATHILDA-FAIL: reason\". The answers are untrusted: Lean \
re-checks them.";

oracleFail[why_String] := "MATHILDA-FAIL: " <> why;

(* Strict text of an evaluated polynomial, no type ascription. *)
oracleText[e_] := With[{x = e}, toLean[HoldForm[x], None, True, {}]];

oracleVarsQ[vars_List] := AllTrue[vars, Head[#] === Symbol &];

Lean`Internal`Oracle["LinearCombination", {a_, b_}, hs_List, vars_List] /;
    oracleVarsQ[vars] && AllTrue[hs, MatchQ[#, {_, _}] &] :=
  Module[{q = Expand[a - b], ps = Expand[#[[1]] - #[[2]]] & /@ hs, cs, txt},
  If[!AllTrue[Append[ps, q], polyQ[#, vars] &],
    Return[oracleFail["the goal and hypotheses must be polynomial in the variables"], Module]];
  cs = Which[q === 0, ConstantArray[0, Length[ps]], ps === {}, $Failed,
    True, cofactors[q, ps, vars]];
  If[cs === $Failed, Return[oracleFail["the goal is not a polynomial combination of the \
hypotheses (no cofactors up to degree 4)"], Module]];
  txt = oracleText /@ cs;
  If[!AllTrue[txt, StringQ], Return[oracleFail["cofactors not translatable to Lean"], Module]];
  "[" <> StringRiffle[txt, ", "] <> "]"];

Lean`Internal`Oracle["SOS", {p_, q_}, gs_List, vars_List] /; oracleVarsQ[vars] :=
  Module[{f = Expand[q - p], gx = Expand /@ gs, r, idx, txt},
  If[!AllTrue[Append[gx, f], polyQ[#, vars] &],
    Return[oracleFail["the goal and hypotheses must be polynomial in the variables"], Module]];
  r = sosTerms[f, gx, vars];
  If[r === $Failed, Return[oracleFail["no sum-of-squares certificate found"], Module]];
  (* 0 for a plain square, i for the multiplier g_i, -1 if mult is neither. *)
  idx[m_] := If[m === 1, 0,
    First[FirstPosition[gx, _?(Expand[# - m] === 0 &), {-1}, {1}]]];
  txt = ({oracleText[#[[1]]], oracleText[#[[2]]], idx[#[[3]]]} &) /@ First[r];
  If[!AllTrue[Flatten[txt[[All, 1 ;; 2]]], StringQ] || MemberQ[txt[[All, 3]], -1],
    Return[oracleFail["squares not translatable to Lean"], Module]];
  "[" <> StringRiffle[("(" <> #[[1]] <> ", " <> #[[2]] <> ", " <>
      ToString[#[[3]]] <> ")" &) /@ txt, ", "] <> "]"];

Lean`Internal`Oracle["Pratt", p_Integer] := Which[
  p === 2, oracleFail["2 is the base of the Pratt scheme (PrattCertifies needs an odd prime)"],
  !PrimeQ[p], oracleFail[ToString[p] <> " is not prime"],
  True, With[{c = pratt[p]}, If[ListQ[c], prattText[c],
    oracleFail["no Lucas witness found for " <> ToString[p]]]]];

Lean`Internal`Oracle["Eval", s_String] := With[{r = Quiet[ToExpression[s]]},
  If[r === $Failed, oracleFail["cannot parse " <> s], ToString[r, LeanForm]]];

Lean`Internal`Oracle[args___] := oracleFail["unsupported query " <> ToString[{args}, InputForm]];

(* ---------------------------------------------------------------------- *)
(* The ProofObject of a Lean (or heuristic) verdict                        *)
(* ---------------------------------------------------------------------- *)

(* The ProofObject association from a generator and a verdict. *)
mkProofObject[g_, v_, engine_, tierShown_, lv_, source_, lt_] := With[{a = <|
  "Kind" -> g["Kind"], "Input" -> g["Input"], "Result" -> g["Result"],
  "Theorem" -> g["Theorem"], "LeanTheorem" -> lt, "ProofScript" -> g["Proof"],
  "LeanSource" -> source, "Certificate" -> g["Certificate"],
  "Verified" -> v["Verified"], "Engine" -> engine, "Tier" -> tierShown,
  "Domain" -> g["Domain"], "Hypotheses" -> g["Hypotheses"],
  "TrustedSteps" -> g["TrustedSteps"], "Axioms" -> v["Axioms"], "Logs" -> v["Logs"],
  "Timing" -> v["Timing"], "Cached" -> TrueQ[v["Cached"]], "LeanVersion" -> lv|>},
  $issued[Hash[a]] = True;
  ProofObject[a]];  (* With: the held slot receives the evaluated Association *)


(* ---------------------------------------------------------------------- *)
(* Encoders: the neutral IR -> Lean text, one rule per obligation "Type"   *)
(* ---------------------------------------------------------------------- *)

(* <|"Statement", "Proof", "Tier", ("Preamble", "UserProof")|>, or $Failed
   after Certify::nocert (the statement does not translate), or
   "Untranslatable" (silent: the IR has a "Fallback"). *)
leanEncode[ir_, o_] := leanEncodeType[ir["Obligation"]["Type"], ir, o];

(* An identity a == b over Int / Rat, or the generic char-0 / ordered-field
   form.  The generic form writes the type variable K for the word Rat: a
   variable printed as K or Rat would be captured by it. *)
leanEncodeType["Identity", ir_, o_] := Module[{ob = ir["Obligation"], gk, lean},
  gk = ob["Generic"];
  If[gk =!= None && AnyTrue[ob["Variables"], MemberQ[{"K", "Rat"}, SymbolName[#]] &],
    Return[nocert[ir["Input"]], Module]];
  lean = toLean[ob["Statement"], ob["Ring"], True, {}, True];
  If[!StringQ[lean], Return[nocert[ir["Input"]], Module]];
  If[gk =!= None, lean = genericStmt[lean, gk]];
  <|"Statement" -> lean, "Proof" -> proofScript[ir["Method"]],
    "UserProof" -> userMethodQ[ir["Method"]], "Tier" -> 0|>];

(* Polynomial identities (one, or a conjunction): intro, then grind on each. *)
leanEncodeType["PolyIdentities", ir_, o_] := Module[{ob = ir["Obligation"], stmt},
  stmt = leanThm[ob["Statement"], ob["Ring"]];
  If[!StringQ[stmt], Return[nocert[ir["Input"]], Module]];
  <|"Statement" -> stmt,
    "Proof" -> "by" <> introLine[ob["Variables"]] <> conjunctionTactic[Length[ob["Pairs"]]],
    "Tier" -> 0|>];
conjunctionTactic[1] := "\n  grind";
conjunctionTactic[n_Integer] := "\n  refine " <> andHoles[n] <> " <;> grind";
andHoles[2] := "And.intro ?_ ?_";
andHoles[n_Integer] := "And.intro ?_ (" <> andHoles[n - 1] <> ")";

(* Ideal membership: have k_i : p_i = 0, have hc : q = Sum c_i * p_i,
   rw [k_i] at hc, grind.  Over every field of characteristic 0 (over Int
   or Q the implication would be weaker than intended: fewer points satisfy
   the hypotheses); the cofactor identity holds in every commutative ring. *)
leanEncodeType["IdealImplication", ir_, o_] := Module[{ob = ir["Obligation"], hs, l, r, cs,
    vars, ty = "Rat", stmt, ptxt, qtxt, ctxt, proof, k},
  hs = ob["Hypotheses"]; {l, r} = ob["Goal"]; cs = ob["Cofactors"]; vars = ob["Variables"];
  stmt = leanThm[ob["Statement"], ty];
  ptxt = leanStmt[inertHeld[iPlus[#[[1]], iTimes[-1, #[[2]]]]], ty] & /@ hs;
  qtxt = leanStmt[inertHeld[iPlus[l, iTimes[-1, r]]], ty];
  ctxt = leanStmt[HoldForm[#], ty] & /@ cs;
  If[!StringQ[stmt] || !AllTrue[Join[ptxt, ctxt, {qtxt}], StringQ],
    Return[nocert[ir["Input"]], Module]];
  k = Length[hs];
  proof = If[MemberQ[ptxt, qtxt], "by" <> introLine[vars] <> "\n  intros\n  grind",
    StringJoin["by", introLine[vars], "\n  intros",
      Table["\n  have mth_k" <> ToString[i] <> " : " <> ptxt[[i]] <> " = 0 := (by grind)",
        {i, k}],
      "\n  have mth_hc : ", qtxt, " = ",
      StringRiffle[Table["(" <> ctxt[[i]] <> ") * (" <> ptxt[[i]] <> ")", {i, k}], " + "],
      " := (by grind)",
      "\n  rw [", StringRiffle[Table["mth_k" <> ToString[i], {i, k}], ", "], "] at mth_hc",
      "\n  grind"]];
  <|"Statement" -> genericStmt[stmt, "CharZero"], "Proof" -> genericText[proof],
    "Tier" -> 0|>];

(* Solve completeness: l - r = factorization, then the equivalence; a text
   that does not translate falls back to the soundness obligation. *)
leanEncodeType["RootsComplete", ir_, o_] := Module[{ob = ir["Obligation"], x, l, r, fact,
    stmt, ftxt, etxt},
  x = ob["Variable"]; {l, r} = ob["Equation"]; fact = ob["Factorization"];
  stmt = leanThm[ob["Statement"], "Rat"];
  ftxt = leanStmt[inertHeld[fact], "Rat"];
  etxt = leanStmt[inertHeld[iEqual[l, r]], "Rat"];
  If[!StringQ[stmt] || !StringQ[ftxt] || !StringQ[etxt], Return["Untranslatable", Module]];
  <|"Statement" -> stmt, "Proof" -> StringJoin["by", introLine[{x}],
      "\n  have mth_e : ", leanStmt[inertHeld[iPlus[l, iTimes[-1, r]]], "Rat"],
      " = ", ftxt, " := (by grind)",
      "\n  have mth_h : (", etxt, ") \[LeftRightArrow] ", ftxt,
      " = 0 := (by rw [\[LeftArrow] mth_e]; grind)",
      "\n  rw [mth_h]\n  grind"], "Tier" -> 0|>];

(* A closed rational fact (Solve soundness, matrices): the kernel evaluates it. *)
leanEncodeType["ClosedFact", ir_, o_] := With[{stmt = leanThm[ir["Obligation"]["Statement"], "Rat"]},
  If[!StringQ[stmt], nocert[ir["Input"]],
    <|"Statement" -> stmt, "Proof" -> "by\n  decide +kernel", "Tier" -> 0|>]];

(* A weighted sum of squares: the hints 0 <= g_i, 0 <= l_k^2 (g), the explicit
   identity target = Sum d_k l_k^2 (g), then grind; every ordered field. *)
leanEncodeType["WeightedSOS", ir_, o_] := Module[{ob = ir["Obligation"], vars, terms, used,
    gnames, stmt, ghints, proof},
  vars = ob["Variables"]; terms = ob["Squares"]; used = ob["Multipliers"];
  gnames = Table[used[[i]] -> "mth_g" <> ToString[i], {i, Length[used]}];
  stmt = leanThm[ob["Statement"], "Rat"];
  If[!StringQ[stmt], Return[nocert[ir["Input"]], Module]];
  ghints = Table["\n  have mth_g" <> ToString[i] <> " : 0 \[LessEqual] " <> ratText[used[[i]]] <>
    " := (by grind)", {i, Length[used]}];
  proof = StringJoin["by", introLine[vars], If[ob["Conditions"] === {}, "", "\n  intros"],
    ghints, sqHints[terms, gnames],
    "\n  have mth_e : ", ratText[ob["Target"]], " = ", If[terms === {}, "0", sqSum[terms]],
    " := (by grind)", "\n  grind"];
  <|"Statement" -> genericStmt[stmt, "Ordered"], "Proof" -> genericText[proof], "Tier" -> 0|>];

(* KKT optimality: the bound for every feasible x (equalities eliminated by
   subst, the multiplier and square hints), then the closed facts. *)
leanEncodeType["KKT", ir_, o_] := Module[{ob = ir["Obligation"], vars, eqHyps, gr, lamFull,
    terms, used, stmt, proof},
  vars = ob["Variables"]; eqHyps = ob["Eliminated"]; gr = ob["Reduced"];
  lamFull = ob["Multipliers"]; terms = ob["Squares"];
  stmt = leanThm[ob["Statement"], "Rat"];
  If[!StringQ[stmt], Return[nocert[ir["Input"]], Module]];
  used = Select[Range[Length[gr]], lamFull[[#]] =!= 0 &];
  proof = StringJoin["by",
    "\n  refine And.intro ?_ ", nestAnd[ob["Facts"]],
    introLine[vars], If[ob["Constraints"] === {}, "", "\n  intros"],
    Table["\n  have mth_v" <> ToString[i] <> " : " <> leanVar[eqHyps[[i, 1]]] <> " = " <>
      ratText[eqHyps[[i, 2]]] <> " := (by grind)\n  subst mth_v" <> ToString[i],
      {i, Length[eqHyps]}],
    Table["\n  have mth_g" <> ToString[j] <> " : 0 \[LessEqual] " <> ratText[gr[[j]]] <>
      " := (by grind)", {j, used}],
    sqHints[terms, {}],
    "\n  have mth_e : ", ratText[ob["Gap"]], " = ",
    StringRiffle[Join[
      ("(" <> ratText[lamFull[[#]]] <> ") * (" <> ratText[gr[[#]]] <> ")") & /@ used,
      If[terms === {}, {}, {sqSum[terms]}]] /. {} -> {"0"}, " + "],
    " := (by grind)", "\n  grind"];
  <|"Statement" -> genericStmt[stmt, "Ordered"], "Proof" -> genericText[proof], "Tier" -> 0|>];

(* Primality: flat Pratt certificates checked by decide +kernel (the
   factorization's odd primes in a preamble of auxiliary theorems), the
   divisors of 2, n < 2, a divisor of a composite. *)
leanEncodeType["Primality", ir_, o_] := Module[{ob = ir["Obligation"], nv, fi, odd, certs, a},
  nv = ob["N"];
  Switch[ob["Form"],
    "Factorization",
      fi = ob["Factorization"]; certs = ob["Certificates"];
      odd = DeleteCases[fi[[All, 1]], 2];
      <|"Preamble" -> StringJoin[Table["theorem @NAME@_pratt_" <> ToString[i] <>
            " : Mathilda.PrattCertifies " <> ToString[odd[[i]]] <> " " <>
            prattText[certs[[i]]] <> " := by decide +kernel\n", {i, Length[odd]}]],
        "Statement" -> StringRiffle[Join[
          {ToString[nv] <> " = " <> StringRiffle[factorPow /@ fi, " * "]},
          ("(\[Exists] c, Mathilda.PrattCertifies " <> ToString[#] <> " c)") & /@ odd],
          " \[And] "],
        "Proof" -> "by\n  exact " <> nestAndTerms[Join[{"by decide +kernel"},
          Table["Exists.intro _ @NAME@_pratt_" <> ToString[i], {i, Length[odd]}]]],
        "Tier" -> 0|>,
    "Small",
      <|"Statement" -> "\[Not] Mathilda.IsPrime " <> ToString[nv],
        "Proof" -> "fun h => absurd h.1 (by decide)", "Tier" -> 0|>,
    "Two",
      <|"Statement" -> "Mathilda.IsPrime 2",
        "Proof" -> "And.intro (by decide) (fun d h => by\n    have hle := Nat.le_of_dvd (by decide) h\n    \
rcases d with _ | _ | _ | d <;>\n      first | exact absurd (Nat.eq_zero_of_zero_dvd h) (by decide) \
| exact Or.inl rfl | exact Or.inr rfl | omega)", "Tier" -> 0|>,
    "Prime",
      <|"Statement" -> "Mathilda.PrattCertifies " <> ToString[nv] <> " " <>
          prattText[First[ob["Certificates"]]],
        "Proof" -> "by decide +kernel", "Tier" -> 0|>,
    "Composite",
      a = ob["Divisor"];
      <|"Statement" -> "\[Not] Mathilda.IsPrime " <> ToString[nv],
        "Proof" -> "fun h => absurd (h.2 " <> ToString[a] <> " (Exists.intro " <>
          ToString[nv/a] <> " (by decide +kernel))) (by decide +kernel)", "Tier" -> 0|>,
    _, nocert[ir["Input"]]]];

(* A closed-form sum: Lean defines t and S by recursion from the lower limit
   a and proves forall m : Nat, S m = closed(m + a) by induction. *)
leanEncodeType["Induction", ir_, o_] := Module[{ob = ir["Obligation"], t, k, av, n, closed,
    tl, cl, ttxt, ctxt, fq},
  t = ob["Term"]; k = ob["Variable"]; av = ob["LowerLimit"]; n = ob["SumVariable"];
  closed = ob["ClosedForm"];
  tl = natForm[t, k, $termBinder, 0];
  cl = natForm[closed, n, $sumBinder, av];
  If[tl === $Failed || cl === $Failed, Return[nocert[ir["Input"]], Module]];
  ttxt = leanStmt[inertHeld[tl], "Rat", $sumPass];
  ctxt = leanStmt[inertHeld[cl], "Rat", $sumPass];
  If[!StringQ[ttxt] || !StringQ[ctxt], Return[nocert[ir["Input"]], Module]];
  (* Fibonacci (v0.6527): the recursive _F first, and in both simp lists;
     every other sum is unchanged. *)
  fq = !FreeQ[{t, closed}, Fibonacci];
  <|"Preamble" -> StringJoin[
      If[fq, "def @NAME@_F : Nat \[RightArrow] Rat\n  | 0 => 0\n  | 1 => 1\n" <>
        "  | n + 2 => @NAME@_F n + @NAME@_F (n + 1)\n", ""],
      "def @NAME@_t (k1 : Nat) : Rat := ", ttxt, "\n",
      "def @NAME@_S : Nat \[RightArrow] Rat\n",
      "  | 0 => @NAME@_t ", ToString[av], "\n",
      "  | m1 + 1 => @NAME@_S m1 + @NAME@_t (m1 + 1", If[av === 0, "", " + " <> ToString[av]],
      ")\n"],
    "Statement" -> "\[ForAll] m1 : Nat, @NAME@_S m1 = " <> ctxt,
    "Proof" -> StringJoin["by\n  intro m1\n  induction m1 with",
      "\n  | zero => simp only [@NAME@_S, @NAME@_t", If[fq, ", @NAME@_F", ""],
      "] <;> (try push_cast) <;> grind",
      "\n  | succ m1 ih => simp only [@NAME@_S, @NAME@_t, ", If[fq, "@NAME@_F, ", ""],
      "ih] <;> (try push_cast) <;> grind"],
    "Tier" -> 0|>];

(* Tier 1: the HasDerivAt walk; Certify::nocert when Tier -> 1 was
   requested, Certify::tier otherwise (no certificate at either tier), and
   none after Certify::vacuous. *)
leanEncodeType["HasDerivAt", ir_, o_] := With[{ob = ir["Obligation"]},
  With[{g = tier1AntiderivGen[ir["Input"], ob["F"], ob["f"], ob["Variable"], o]},
    Which[AssociationQ[g], g, g === "vacuous", $Failed, True,
      If[o["Tier"] === 1, Message[Certify::nocert, ir["Input"]],
        Message[Certify::tier, ir["Input"]]];
      $Failed]]];

leanEncodeType[_, ir_, o_] := nocert[ir["Input"]];

(* ---------------------------------------------------------------------- *)
(* The "Lean" backend: runs and presentation                               *)
(* ---------------------------------------------------------------------- *)

(* One theorem chunk per item, named after its position in the call. *)
leanChunks[items_List, o_] := With[{hb = heartbeats[o["TimeConstraint"]]},
  theoremChunk["mathilda_cert_" <> ToString[#["Index"]], #["G"]["Statement"],
    #["G"]["Proof"], hb, Lookup[#["G"], "Preamble", ""]] & /@ items];

(* The item's generator association with its auxiliary names (@NAME@) under
   the theorem's actual name, as the ProofObject shows them (the cache keys
   use the generator as it is). *)
leanNamed[item_] := With[{g = item["G"], nm = "mathilda_cert_" <> ToString[item["Index"]]},
  Join[g, AssociationMap[StringReplace[g[#], "@NAME@" -> nm] &,
    Select[{"Statement", "Proof", "Preamble"}, StringQ[g[#]] &]]]];

(* Per item, the source of the file its tier's items form. *)
leanBatchSources[items_List, o_] := Module[{chunks = leanChunks[items, o], idx,
    src = ConstantArray[None, Length[items]]},
  Do[
    idx = Select[Range[Length[items]], items[[#]]["G"]["TierUsed"] === tier &];
    If[idx =!= {}, With[{s = First[assembleFile[chunks[[idx]], tier]]},
      Do[src[[j]] = s, {j, idx}]]],
    {tier, {0, 1}}];
  src];

(* One Lean run per tier (cached verdicts are reused; items with a user
   Method script get a file of their own, so no chunk can influence another
   one's verdict).  Verify -> False: no run, Missing["NotChecked"]. *)
leanRun[items_List, o_] := Module[{gs = #["G"] & /@ items, chunks, out, hb, exe, proj, idx,
    keys, todo, fresh, v, lv},
  hb = heartbeats[o["TimeConstraint"]];
  chunks = leanChunks[items, o];
  out = ConstantArray[None, Length[items]];
  proj = leanProject[o["LeanProject"]];
  exe = leanExe[o["LeanExecutable"]];
  Do[
    idx = Select[Range[Length[items]], gs[[#]]["TierUsed"] === tier &];
    Which[
      idx === {}, Null,
      !TrueQ[o["Verify"]],
        With[{src = First[assembleFile[chunks[[idx]], tier]]},
          Do[out[[j]] = <|"Verified" -> Missing["NotChecked"], "Axioms" -> {}, "Logs" -> {},
            "Timing" -> 0., "Source" -> src, "LeanVersion" -> Missing["NotChecked"]|>,
            {j, idx}]],
      True,
        keys = AssociationThread[idx -> (cacheKey[gs[[#]], hb, exe, proj] & /@ idx)];
        todo = Select[idx, !KeyExistsQ[$certCache, keys[#]] &];
        fresh = <||>;
        If[todo =!= {},
          With[{shared = Select[todo, !TrueQ[gs[[#]]["UserProof"]] &],
              own = Select[todo, TrueQ[gs[[#]]["UserProof"]] &]},
            If[shared =!= {},
              fresh = AssociationThread[shared -> runBatch[chunks[[shared]], tier, exe, proj]]];
            Do[fresh = Append[fresh, j -> First[runBatch[chunks[[{j}]], tier, exe, proj]]],
              {j, own}]];
          Do[If[BooleanQ[fresh[j]["Verified"]],
              $certCache = Append[$certCache, keys[j] -> fresh[j]]], {j, todo}]];
        Do[
          v = If[KeyExistsQ[fresh, j], Append[fresh[j], "Cached" -> False],
            Append[$certCache[keys[j]], "Cached" -> True]];
          lv = If[MissingQ[v["Verified"]], Missing["NotAvailable"], leanVersionOf[exe]];
          If[v["Verified"] === False, Message[Certify::rejected, gs[[j]]["Input"]]];
          out[[j]] = Append[v, "LeanVersion" -> lv],
          {j, idx}]],
    {tier, {0, 1}}];
  out];

leanPresent[item_, v_, o_] := With[{g = leanNamed[item]},
  mkProofObject[g, v, "Lean", g["TierUsed"], v["LeanVersion"], v["Source"], g["Statement"]]];

(* A refutation (Refute.m): the theorem not (S), S the statement text as the
   footer states it (generic binders included, so the theorem itself has
   none), proved from the witness:
     fun h => absurd (h [(K := Rat)] w1 ... [(by decide +kernel) per
       antecedent]) (by decide +kernel)
   a closed statement by decide +kernel, a sum by h m0 (m0 = n0 - a, the
   preamble's definitions reused).  The result is a generator association
   for leanRun. *)
leanRefuteEncode[rir_, o_] := Module[{p = rir["Problem"], w = rir["Witness"], ty, gen, proof},
  ty = If[p["DomainForm"] === "Int", "Int", "Rat"];
  gen = MemberQ[{"CharZeroField", "OrderedField"}, p["DomainForm"]];
  proof = Switch[p["Form"],
    "Closed", "by\n  decide +kernel",
    "Induction", "fun h => absurd (h " <> ToString[w[[1, 2]] - p["Shift"]] <>
      ") (by decide +kernel)",
    _, StringJoin["fun h => absurd (h", If[gen, " (K := Rat)", ""],
      (" (" <> refValueText[#, ty] <> " : " <> ty <> ")") & /@ w[[All, 2]],
      ConstantArray[" (by decide +kernel)", p["Antecedents"]], ") (by decide +kernel)"]];
  <|"Statement" -> "\[Not] (" <> rir["LeanStatement"] <> ")", "Proof" -> proof,
    "Preamble" -> rir["Preamble"], "Tier" -> 0, "TierUsed" -> 0, "UserProof" -> False,
    "Input" -> rir["Input"]|>];
refValueText[v_, ty_] := With[{s = toLean[HoldForm[v], ty, True, {}]}, If[StringQ[s], s, "0"]];

registerBackend["Lean", <|"Name" -> "Lean", "Trust" -> "Kernel", "Encode" -> leanEncode,
  "Run" -> leanRun, "Present" -> leanPresent, "Refute" -> leanRefuteEncode,
  "Label" -> ("Lean Tier " <> ToString[#["Tier"]] &),
  "Options" -> {"Tier", "Method", "TimeConstraint", "LeanExecutable", "LeanProject"}|>];

End[];

EndPackage[];
