(* ====================================================================== *)
(* Lean.m -- Lean 4 certification: Certify, ProofObject, LeanVerify.       *)
(*                                                                          *)
(* Loaded at startup by core_init (after FullSimplify.m and                 *)
(* SeparablePolynomialDecomposition.m).  Defines rules only: no Lean process *)
(* is started and no file is read at load time.                             *)
(*                                                                          *)
(* Trust model.  Mathilda is an untrusted oracle and the Lean kernel is the *)
(* trusted checker.  Certify computes an answer and a certificate, emits a  *)
(* self-contained Lean 4 file, and reports "Verified" -> True only when     *)
(* Lean accepted the theorem with no sorry and no axiom outside             *)
(* {propext, Classical.choice, Quot.sound}.  "TrustedSteps" lists, per      *)
(* certificate, the facts Lean did NOT check (e.g. the Together normal form *)
(* used to clear denominators).                                             *)
(*                                                                          *)
(* The C side (src/print_lean.c, src/lean_run.c) supplies Lean`Internal`ToLean (the  *)
(* strict LeanForm translator), Lean`Internal`RunProcess and Lean`Internal`VerifySource (the  *)
(* system()-based runner and NDJSON parser).                                *)
(* ====================================================================== *)

BeginPackage["Lean`"];

Certify::usage =
  "Certify[comp] evaluates the computation comp (Expand, Factor, FactorTerms, Collect, \
Together, Apart, Cancel, Simplify, FullSimplify, Integrate over Q(x), \
PolynomialExtendedGCD, PolynomialReduce, Solve, Inverse, Det, LinearSolve, NullSpace, \
Minimize, Maximize, FactorInteger, PrimeQ, Sum) and certifies its result with the Lean 4 \
kernel; Certify[a == b] certifies the claim a == b, Certify[Implies[p1 == q1 && ..., \
p == q]] an ideal membership and Certify[ForAll[vars, p >= q]] a nonnegativity (sum of \
squares). Tier -> 1 gives Mathlib HasDerivAt antiderivatives. Certify[{c1, c2, ...}] \
certifies a batch in one Lean run. The result is a ProofObject whose \"Verified\" property \
is True only when Lean accepted the theorem.";

ProofObject::usage =
  "ProofObject[<|...|>] is a Lean 4 certificate produced by Certify; obj[\"prop\"] gives a \
property and obj[\"Properties\"] lists them; obj[\"Verified\"] is \
Missing[\"Unauthenticated\"] for an object Certify did not produce in this session. \
ProofObject[comp, opts] (comp not an Association) is the same as Certify[comp, opts].";

LeanVerify::usage =
  "LeanVerify[src] checks the Lean 4 source string src with the Lean executable and gives \
<|\"Success\", \"ExitCode\", \"Messages\", \"Axioms\", \"AxiomReports\", \"Timing\", \
\"Command\"|>. Success requires exit code 0, no error, no sorry, and only the standard axioms \
for every declaration (examples are named so they are checked too). Only sources with no \
user-defined code are run: no imports (core Lean) or only import Mathlib with a project, and \
only theorem / example / section / variable / universe / numeric \
set_option commands; anything else (def, instance, macro, attribute, #eval, native_decide, \
sorry, ...) is refused unrun.";

$LeanExecutable::usage =
  "$LeanExecutable is the Lean 4 executable used by Certify and LeanVerify: Automatic \
(MATHILDA_LEAN, else the newest ~/.elan toolchain, else lean on PATH) or a path.";

$LeanProject::usage =
  "$LeanProject is None or the directory of a Lake project with Mathlib built, used for \
Tier -> 1 certificates. It defaults to the environment variable MATHILDA_LEAN_PROJECT.";

$LeanVersion::usage =
  "$LeanVersion gives the version string of the Lean executable, or \
Missing[\"NotAvailable\"] when none is found.";

Engine::usage =
  "Engine is an option of Certify and ProofObject: \"Lean\" (the Lean kernel judges) or \
\"Mathilda\" (a heuristic PossibleZeroQ / Simplify check, for contrast; its \"Verified\" is \
Missing[\"Heuristic\", True | False], never True).";

Tier::usage =
  "Tier is an option of Certify: 0 (core Lean), 1 (Mathlib, needs $LeanProject) or \
Automatic.";

Domain::usage =
  "Domain is an option of Certify for identity claims: Automatic (Mathematica's default, \
complex variables), Complexes (every field of characteristic 0), Reals (every ordered \
field), Rationals (Rat) or Integers (Int); any other value is refused (Certify::domain).";

Verify::usage =
  "Verify is an option of Certify: False generates the Lean source without running Lean \
(\"Verified\" -> Missing[\"NotChecked\"]).";

LeanExecutable::usage =
  "LeanExecutable is an option of Certify and LeanVerify that overrides $LeanExecutable.";

LeanProject::usage =
  "LeanProject is an option of Certify and LeanVerify that overrides $LeanProject.";

Oracle::usage =
  "Oracle[kind, args...] is the server side of the Lean 4 tactics in the lean/ package \
(mathilda_linear_combination, mathilda_sos, mathilda_prime, #mathilda). It answers with \
plain data as a String (numerals, v1, v2, ..., + - * / ^, lists and tuples), which the Lean \
side parses itself and never elaborates as Lean syntax: Oracle[\"LinearCombination\", {a, b}, {{l1, r1}, ...}, vars] gives \
cofactors c with a - b == Sum c_i (l_i - r_i); Oracle[\"SOS\", {p, q}, {g1, ...}, vars] gives \
(d, l, i) triples with q - p == Sum d l^2 g_i (g_0 = 1); Oracle[\"Pratt\", p] gives a flat \
Pratt certificate; Oracle[\"Eval\", \"input\"] gives the LeanForm of the evaluated input. \
A failure gives \"MATHILDA-FAIL: reason\". The answers are untrusted: Lean re-checks them.";

KeepFiles::usage =
  "KeepFiles is an option of LeanVerify: True keeps the temporary .lean file \
(its path is returned as \"SourceFile\").";

Begin["`Private`"];

(* ---------------------------------------------------------------------- *)
(* Messages                                                                *)
(* ---------------------------------------------------------------------- *)

LeanVerify::nolean = "no Lean 4 executable found; set $LeanExecutable or MATHILDA_LEAN";
LeanVerify::badpath = "path `1` cannot be passed safely to the shell on this platform";
LeanVerify::fail = "Lean process could not be started: `1`";
Certify::nocert = "no Lean certificate is available for `1`";
Certify::engine = "unknown Engine `1`; use \"Lean\" or \"Mathilda\"";
Certify::tier = "`1` needs Tier 1 (Mathlib); set $LeanProject";
Certify::nomathlib = "Tier -> 1 requested but $LeanProject is not set";
Certify::rejected = "Lean rejected the certificate for `1`; see obj[\"Logs\"]";
Certify::method = "Method `1` is not accepted: `2`; give a plain tactic script (see ?Certify)";
LeanVerify::unsafe = "the source uses `1`, which can run code at elaboration time, bypass or forge the kernel check; not run";
LeanVerify::import = "the source has `1`; LeanVerify accepts no imports (core Lean) or, with a LeanProject, only import Mathlib / import Mathlib.<Path>; not run";
LeanVerify::cmd = "the source uses the command `1`; LeanVerify accepts only theorem, example, section, end, variable, universe and numeric / Boolean set_option commands (lemma with Mathlib), each starting at column 0; not run";
LeanVerify::lex = "line `2`: the source contains `1`, which is outside LeanVerify's minimal lexical surface (no strings, comments, character literals, #, $, %, backquotes, backslashes, @[, French quotes or unlisted Unicode); not run";
Certify::unsafe = "the generated Lean source failed the LeanVerify source check (`1`); not run";
Certify::eval = "the claim `1` applies `2`, which Mathilda would have to evaluate before Lean sees the claim; certify the computation itself (e.g. Certify[Simplify[...]]) or state the claim with arithmetic only";
Certify::value = "`1` has `2`; Clear[`1`], or certify the value-substituted claim explicitly (write it out, or Certify[Evaluate[claim]]); a claim is certified as written, never with values substituted";
Certify::context = "`1` is in the context `2`, which is reserved for Certify's own generator symbols; use a symbol of your own";
Certify::domain = "Domain -> `1` is not supported; use Automatic, Complexes, Reals, Rationals or Integers";
Certify::element = "the assumption `1` cannot be stated for Domain -> `2`; only Element[vars, Reals] / Element[vars, Complexes] (and Element[vars, d] implied by the domain d) are accepted";
Certify::optdefault = "the default Assumptions -> `1` was set with SetOptions[Certify, ...], which evaluated it before Certify could check it; pass Assumptions in the call instead, or reset it with SetOptions[Certify, Assumptions -> {}]";
Certify::antideriv = "`1` is certified as `2`: the claimed F is an antiderivative of f; the equation itself holds only up to an additive locally constant term (see \"TrustedSteps\")";
Certify::realvars = "the claim `1` has order relations, but `2` is not constrained to be real; add Element[`2`, Reals] (or an order relation on it) to the assumptions";
Certify::optignored = "`1` is ignored for a `2` certificate (only identity claims read Domain, Assumptions and Method)";

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

Options[LeanVerify] = {LeanExecutable -> Automatic, LeanProject -> Automatic,
  KeepFiles -> False};

LeanVerify[src_String, opts:OptionsPattern[]] := Module[{exe, proj, full, lines, r, chk,
    ok},
  exe = leanExe[OptionValue[LeanExecutable]];
  proj = leanProject[OptionValue[LeanProject]];
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
  r = Lean`Internal`VerifySource[full, exe, proj, TrueQ[OptionValue[KeepFiles]],
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

leanType["Int"] = Integers; leanType["Rat"] = Rationals; leanType["Real"] = Reals;

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

(* ---------------------------------------------------------------------- *)
(* Certificate generators (no Lean call)                                   *)
(* ---------------------------------------------------------------------- *)

(* Domain selection.  Int is chosen automatically only for a hypothesis-free
   polynomial identity with integer coefficients: such an identity holding
   for all integers holds as a polynomial identity, hence in every
   commutative ring, so Int loses nothing.  With hypotheses (Assumptions,
   implications, inequalities) quantifying over Int would weaken the claim --
   it could even be vacuous (forall x : Int, 2 x = 1 -> ...) -- so those
   statements are stated over Rat unless Domain -> Integers is given. *)
autoDomain[stmt_] := If[FreeQ[stmt, _Rational] && FreeQ[stmt, Power[_, _?Negative]],
  "Int", "Rat"];

(* The Domain option (sixth review, v0.6522).  Each accepted value has one
   meaning, applied by identityGen (the only generator that reads it):
     Automatic  Mathematica's default (variables complex unless an order
                relation or Element[v, Reals] makes them real): a
                hypothesis-free identity over Int (integer polynomial) or Rat
                (equivalent: a polynomial identity with rational
                coefficients holding on all of Z or Q holds in Q[x], hence in
                every field of characteristic 0); with hypotheses the generic
                char-0 / ordered-field form, as for Complexes.
     Complexes  every field of characteristic 0 (C among them), hypotheses or
                not; the ordered-field form when every variable is declared
                real; order relations on a variable that is not real are
                refused (Certify::realvars).
     Reals      every ordered field (R among them): all variables real by
                declaration.
     Rationals  Rat.      Integers  Int.
   Anything else is refused up front (Certify::domain). *)
$domains = {Automatic, Complexes, Reals, Rationals, Integers};

(* Generic domains.  Mathematica's variables are complex unless an inequality
   makes them real, so a statement with hypotheses is proved for EVERY field
   of characteristic 0 (it then holds in Q, R and C alike), or, when it
   involves an order, for EVERY ordered field (R among them).  A hypothesis
   that only Q cannot satisfy (x^2 == 2, x^2 + 4 == 0) no longer makes the
   statement vacuous: forall x, x^2 + 4 = 0 -> x = 7 is false in C, so Lean
   cannot prove it.  The statement is first translated over Rat and then
   generalised: the word Rat becomes the type variable K, bound by
   $charZeroBinders / $orderedBinders (core Lean's grind classes).  The
   theorem takes them as parameters; the footer states the full
   forall {K : Type} [...], ... form. *)
$charZeroBinders = "{K : Type} [Lean.Grind.Field K] [Lean.Grind.IsCharP K 0]";
$orderedBinders = "{K : Type} [Lean.Grind.Field K] [LE K] [LT K] [Std.LawfulOrderLT K] \
[Std.IsLinearOrder K] [Lean.Grind.OrderedRing K]";
$genericDomain["CharZero"] = Complexes; $genericDomain["Ordered"] = Reals;
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

(* An order relation anywhere in e makes the statement real. *)
orderedQ[e_] := !FreeQ[e, Less | Greater | LessEqual | GreaterEqual | iLess | iGreater |
  iLessEqual | iGreaterEqual];

(* The variables of relations (x > 0, a == b, ...) and expressions. *)
relationVariables[e_] := Variables[Flatten[{e} //. {HoldPattern[Inequality[a___]] :> Take[{a}, {1, -1, 2}],
    (Equal | Unequal | Less | LessEqual |
    Greater | GreaterEqual | And | Or | Not | Implies)[a___] :> {a}}]];

polyQ[e_, vars_] := PolynomialQ[e, vars] && FreeQ[e, Power[_, _?Negative]];

(* {numerator, denominator} of e read off its structure, with NO cancellation
   (Together would cancel common factors before Lean sees anything, turning
   (x^2 - 1)/(x - 1) == x + 1 into 1 + x = 1 + x).  Sums are put over the
   product of their denominators. *)
rawFraction[e_Plus] := Fold[{#1[[1]] #2[[2]] + #2[[1]] #1[[2]], #1[[2]] #2[[2]]} &,
  rawFraction /@ List @@ e];
rawFraction[e_Times] := Fold[{#1[[1]] #2[[1]], #1[[2]] #2[[2]]} &, {1, 1},
  rawFraction /@ List @@ e];
rawFraction[Power[b_, k_Integer]] := With[{f = rawFraction[b]},
  If[k >= 0, f^k, Reverse[f]^(-k)]];
rawFraction[e_] := {e, 1};

(* The same, read off a HELD side (before Mathilda's arithmetic could cancel
   x/x to 1 or x^2/x to x): each numerator and denominator is a polynomial
   expression, whose evaluation is ring arithmetic only. *)
SetAttributes[{rawH, rawHl}, HoldAllComplete];
rawH[Plus[ts__]] := Fold[{#1[[1]] #2[[2]] + #2[[1]] #1[[2]], #1[[2]] #2[[2]]} &,
  {rawHl[ts]}];
rawH[Times[fs__]] := Fold[{#1[[1]] #2[[1]], #1[[2]] #2[[2]]} &, {1, 1}, {rawHl[fs]}];
rawH[Power[b_, k_Integer]] := With[{f = rawH[b]}, If[k >= 0, f^k, Reverse[f]^(-k)]];
rawH[Subtract[a_, b_]] := rawH[Plus[a, Times[-1, b]]];
rawH[Divide[a_, b_]] := rawH[Times[a, Power[b, -1]]];
rawH[Minus[a_]] := rawH[Times[-1, a]];
rawH[e_] := rawFraction[e];
rawHl[] := Sequence[];
rawHl[t_, r___] := Sequence[rawH[t], rawHl[r]];

(* The certified statement for the identity a == b over `vars`:
   {held statement, hypotheses, form, extra certificate data, extra trusted steps}.
   A rational identity is cross-multiplied with the unreduced numerators and
   denominators of the two sides, so Lean checks the actual polynomial
   identity n_a d_b = n_b d_a. *)
identityStatement[a_, b_, vars_, raw_:Automatic] := Module[{na, da, nb, db},
  {{na, da}, {nb, db}} = If[raw === Automatic, {rawFraction[a], rawFraction[b]}, raw];
  (* A nonzero numeric denominator is no side condition: fold it back. *)
  If[NumberQ[da] && da =!= 0, na = na/da; da = 1];
  If[NumberQ[db] && db =!= 0, nb = nb/db; db = 1];
  If[da === 1 && db === 1 && polyQ[a, vars] && polyQ[b, vars],
    Return[{HoldForm[a == b], {}, "Polynomial", <||>, {}}, Module]];
  If[!AllTrue[{na, da, nb, db}, polyQ[#, vars] &] || Expand[da] === 0 || Expand[db] === 0,
    Return[$Failed, Module]];
  With[{l = If[db === 1, na, na db], r = If[da === 1, nb, nb da]},
    {HoldForm[l == r],
     Join[If[da === 1, {}, With[{d = da}, {HoldForm[d != 0]}]],
          If[db === 1, {}, With[{d = db}, {HoldForm[d != 0]}]]],
     "ClearedDenominators",
     <|"Numerators" -> {na, nb}, "Denominators" -> {da, db}|>,
     {"field algebra: the cross-multiplied polynomial identity (numerators and \
denominators of each side as written, no cancellation) gives the rational identity \
wherever every denominator is nonzero (listed in \"Hypotheses\")"}}]];

(* Real variables (finding B, v0.6520).  The ordered-field form
   (forall ordered fields K, forall x y : K, ...) is faithful to
   Mathematica's semantics only when every bound variable is real: in
   Mathematica an order relation makes ITS variables real, not the others
   (x^2 + 4 == 0 && y > 0 has the complex solution x = 2 I).  A variable is
   real when it occurs in an order relation of the assumptions or is
   declared with Element[v, Reals] / Element[{v, ...}, Reals]. *)
$orderHeads = Less | Greater | LessEqual | GreaterEqual;
orderVariables[assum_List] := Union[Variables[Join[
  Cases[assum, ($orderHeads)[xs__] :> {xs}, {0, Infinity}],
  Cases[assum, HoldPattern[Inequality[xs___]] :> Take[{xs}, {1, -1, 2}], {0, Infinity}]]]];
(* {other assumptions, variables declared real} of a list of assumptions, or
   $Failed after Certify::element.  The same rules hold for every Domain:
   Element[v, Reals] declares v real (a generic domain then uses the
   ordered-field form; over Int / Rat it is implied and dropped);
   Element[v, Complexes] is implied by every domain and dropped (dropping a
   hypothesis can only strengthen the statement); Element[v, Rationals] is
   implied over Rat and Int, Element[v, Integers] over Int; any other
   Element cannot be stated and is refused. *)
splitRealDecls[assum_List, dom_:Automatic] := Module[{implied, bad},
  implied = Switch[dom, Integers, Integers | Rationals | Reals | Complexes,
    Rationals, Rationals | Reals | Complexes, _, Complexes];
  bad = Select[assum, MatchQ[#, Element[_, _]] && !MatchQ[#, Element[_, Reals]] &&
      !MatchQ[#, Element[_, implied]] &];
  If[bad =!= {}, Message[Certify::element, First[bad], dom]; Return[$Failed, Module]];
  {DeleteCases[assum, Element[_, Reals | implied]],
   If[MemberQ[{Integers, Rationals}, dom], {},
     Union[Flatten[Cases[assum, Element[v_, Reals] :> {v}]]]]}];

(* Generator for an identity a == b (rows 1-2).  Returns the generator
   association, or $Failed after a message. *)
identityGen[kind_, input_HoldForm, result_, a_, b_, trusted_List, o_, raw_:Automatic] :=
  Module[{vars, st, ty, held, assum, lean, tier, gk, reals, ordered, notReal, dom, sp, auto},
  dom = o["Domain"];
  sp = splitRealDecls[DeleteCases[Flatten[{o["Assumptions"]} /. And -> List], True], dom];
  If[sp === $Failed, Return[$Failed, Module]];
  {assum, reals} = sp;
  If[!AllTrue[reals, Head[#] === Symbol &],
    Message[Certify::nocert, input]; Return[$Failed, Module]];
  vars = Sort[Union[Variables[{a, b}], If[ListQ[raw], Variables[raw], {}],
    relationVariables[assum], reals]];
  If[!AllTrue[vars, Head[#] === Symbol &],
    Message[Certify::nocert, input]; Return[$Failed, Module]];
  st = identityStatement[a, b, vars, raw];
  If[st === $Failed, Message[Certify::nocert, input]; Return[$Failed, Module]];
  tier = 0;
  (* The statement's domain ($domains above).  Mathematica's semantics
     (Automatic with hypotheses, Complexes): every char-0 field (variables
     complex), or every ordered field when every variable is real; order
     relations with a variable that is not real: no certificate. *)
  ordered = orderedQ[{assum, st[[2]]}];
  notReal = Complement[vars, orderVariables[assum], reals];
  auto := Which[
    ordered && notReal =!= {},
      Message[Certify::realvars, input, If[Length[notReal] === 1, First[notReal], notReal]];
      $Failed,
    ordered || (reals =!= {} && notReal === {}), "Ordered",
    True, "CharZero"];
  gk = Switch[dom,
    Automatic, If[assum === {}, None, auto],
    Complexes, auto,
    Reals, "Ordered",
    _, None];
  If[gk === $Failed, Return[$Failed, Module]];
  (* The generic form writes the type variable K for the word Rat: a
     variable printed as K or Rat would be captured by it. *)
  If[gk =!= None && AnyTrue[vars, MemberQ[{"K", "Rat"}, SymbolName[#]] &],
    Message[Certify::nocert, input]; Return[$Failed, Module]];
  ty = Switch[dom, Integers, "Int", Automatic, If[gk === None, autoDomain[First[st]], "Rat"],
    _, "Rat"];
  (* Bind the variables (and assumptions) around the held statement without
     ever evaluating it: Equal of two identical sides would become True. *)
  held = With[{vs = vars, h = If[Length[assum] === 1, First[assum], And @@ assum]},
    Which[
      vs === {} && assum === {}, First[st],
      vs === {}, First[st] /. HoldForm[s_] :> HoldForm[Implies[h, s]],
      assum === {}, First[st] /. HoldForm[s_] :> HoldForm[ForAll[vs, s]],
      True, First[st] /. HoldForm[s_] :> HoldForm[ForAll[vs, h, s]]]];
  lean = toLean[held, ty, tier === 0, {}, True];
  If[!StringQ[lean], Message[Certify::nocert, input]; Return[$Failed, Module]];
  If[gk =!= None, lean = genericStmt[lean, gk]];
  <|"Kind" -> kind, "Input" -> input, "Result" -> result, "Theorem" -> held,
    "Statement" -> lean, "Proof" -> proofScript[o["Method"]],
    "UserProof" -> userMethodQ[o["Method"]],
    "Certificate" -> Join[<|"Form" -> st[[3]], "Variables" -> vars|>, st[[4]]],
    "TrustedSteps" -> Join[trusted, st[[5]]], "Hypotheses" -> st[[2]],
    "Tier" -> tier, "Domain" -> If[gk =!= None, $genericDomain[gk], leanType[ty]],
    "Claim" -> {a, b}|>];

$identityHeads = {Expand, Collect, Together, Apart, Cancel, Simplify, FullSimplify};
$factorHeads = {Factor, FactorTerms};

(* ---- Inputs are translated as written (finding C, v0.6520) -----------------
   Certify used to evaluate both sides of a claim before translating them, so
   Simplify[(x^2 - 1)/(x - 1)] == x + 1 became the tautology 1 + x = 1 + x
   (likewise Integrate, D, Refine, Limit, Mod, Exp[Log[x]], ...), which Lean
   "verified" with no trusted step.  Now every input is walked HELD
   (rewriteInput) before any generator sees it:
     - numbers, symbols, Plus / Times / Subtract / Divide / Minus, integer
       powers and the relations / connectives / List of a claim are kept
       (Mathilda's evaluation of them is arithmetic normalisation, which Lean
       re-checks: the statement is a polynomial identity);
     - Expand, Together, Cancel, Apart, Factor, FactorTerms and Collect are
       value-preserving on rational functions, so h[e] is replaced by e: Lean
       then checks the user's own expression (a Factor / FactorTerms adds the
       trusted step "irreducibility of the factors");
     - D[F, x] with F a polynomial is replaced by Mathilda's derivative, and
       "Mathilda's differentiation" is a trusted step;
     - an indefinite Integrate[f, x] as a whole side of a claim goes to the
       antiderivative generator with the CLAIMED antiderivative (Lean checks
       D[F] == f through the cleared chain-rule identity); a Sum[t, {k, a, n}]
       side goes to the closed-form sum generator with the claimed closed form;
     - the outermost head of a computation (Certify[Simplify[e]]) is kept;
     - anything else (Simplify, Refine, Limit, Mod, Exp, Log, Sqrt, a definite
       Integrate, a user function, ...) inside the input: Certify::eval. *)
(* HoldAllComplete (finding 4, v0.6521): with HoldAll, an Evaluate[...]
   inside the claim was evaluated as soon as the walker received it, so
   Evaluate[Simplify[(x^2 - 1)/(x - 1)]] == x + 1 became a tautology.  Under
   HoldAllComplete nothing in the claim is evaluated by the walker, and
   Evaluate itself is refused (Certify::eval) like any other head.  A symbol
   with a value (x = 5) is refused with Certify::value (finding 5): the claim
   is certified as written, never with values substituted. *)
SetAttributes[{cw, cwl, cwClaim, hk, csym, valuedQ, symbolProblem, cvalue}, HoldAllComplete];
(* Why the (held) symbol s cannot stand for itself in a claim or in the
   Assumptions, or None (sixth review, v0.6522; before, Protected symbols
   were exempt and only OwnValues were looked at, so Protect[e] or an
   UpValue t /: t^2 = -1 turned a claim into a tautology).
     - Context Lean`Private` / Lean`Internal`: Certify's own generator
       symbols (surd variables, prelude placeholders, sum binders); a claim
       naming one could capture them -> {"context", ctx}.
     - System` symbols are the built-in heads and constants the walker and
       the strict translator already understand (Plus, Power, Pi, I,
       Reals, True, ...); they are accepted unless one has an OwnValue or an
       UpValue a user gave it (I's own OwnValue I :> I is not a value).
     - Every other symbol must be a plain variable: no OwnValues,
       DownValues, UpValues, SubValues, NValues or DefaultValues, and no
       attribute that changes how it evaluates or differentiates
       (Constant, NumericFunction, Flat, Orderless, Listable, Hold*, ...);
       Protected, ReadProtected, Locked and Temporary are harmless. *)
$systemValueAllowed = {"I"};  (* names: I itself evaluates to Complex[0, 1] *)
$inertAttributes = {Protected, ReadProtected, Locked, Temporary};
symbolProblem[s_Symbol] := Module[{ctx = Context[s], k},
  Which[
    StringStartsQ[ctx, "Lean`Private`" | "Lean`Internal`"], {"context", ctx},
    ctx === "System`",
      Which[
        MemberQ[$systemValueAllowed, SymbolName[Unevaluated[s]]], None,
        OwnValues[s] =!= {}, {"value", "a value (" <> ToString[OwnValues[s][[1, 2]], InputForm] <> ")"},
        UpValues[s] =!= {}, {"value", "UpValues (" <> ToString[UpValues[s], InputForm] <> ")"},
        True, None],
    True,
      k = SelectFirst[{{"a value", OwnValues[s]}, {"DownValues", DownValues[s]},
          {"UpValues", UpValues[s]}, {"SubValues", SubValues[s]}, {"NValues", NValues[s]},
          (* DefaultValues only where it is a HoldAll built-in: never evaluate s. *)
          {"DefaultValues", If[MemberQ[Attributes[DefaultValues], HoldAll],
            With[{d = Quiet[DefaultValues[s]]}, If[ListQ[d], d, {}]], {}]},
          {"the attributes", Complement[Attributes[s], $inertAttributes]}},
        #[[2]] =!= {} &, None];
      Which[
        k === None, None,
        k[[1]] === "a value", {"value", "a value (" <> ToString[k[[2, 1, 2]], InputForm] <> ")"},
        True, {"value", k[[1]] <> " (" <> ToString[k[[2]], InputForm] <> ")"}]]];
symbolProblem[_] := None;
valuedQ[s_] := symbolProblem[s] =!= None;
symbolMessage[s_HoldForm, {"context", ctx_}] := Message[Certify::context, s, ctx];
symbolMessage[s_HoldForm, {_, what_}] := Message[Certify::value, s, what];
$cwKeep = {Plus, Times, Subtract, Divide, Minus, Equal, Unequal, Less, LessEqual, Greater,
  GreaterEqual, And, Implies, Inequality, ForAll, List, Element};
$cwComputations = {Expand, Collect, Together, Apart, Cancel, Simplify, FullSimplify,
  Factor, FactorTerms, Integrate, PolynomialExtendedGCD, PolynomialReduce, Solve, Inverse,
  Det, LinearSolve, NullSpace, Minimize, Maximize, FactorInteger, PrimeQ, Sum};
$cwAnti = False;    (* inside a claimed antiderivative: surds, Log, ArcTan, ArcTanh *)
$cwAnyPow = False;  (* inside a sum: c^k with a symbolic exponent *)
$cwAssum = False;   (* walking the Assumptions option *)
$cwNotes = {};      (* heads stripped / trusted steps added during one walk *)
$cwTrusted = {};

cw[e_Integer] := e;
cw[e_Rational] := e;
cw[e_Real] := e;
cw[e_String] := e;
cw[s_Symbol] /; valuedQ[s] := cvalue[s, symbolProblem[s]];
cw[s_Symbol] := hk[s];
(* A variable slot (the x of Integrate[f, x], the k of a Sum) is checked
   like any other symbol of the claim. *)
csym[s_Symbol] := If[valuedQ[s], cvalue[s, symbolProblem[s]], hk[s]];
cw[Evaluate[___]] := cfail[Evaluate];
cw[Power[b_, n_Integer]] := inert[Power, cw[b], n];
cw[Power[b_, e_]] /; $cwAnti || $cwAnyPow := inert[Power, cw[b], cw[e]];
cw[(h:Sqrt | Log | ArcTan | ArcTanh)[u_]] /; $cwAnti := inert[h, cw[u]];
cw[(h:Expand | Together | Cancel | Factor | FactorTerms)[e_]] :=
  ($cwNotes = Append[$cwNotes, h]; cw[e]);
cw[(h:Apart | Collect)[e_, v_]] /; MatchQ[Unevaluated[v], _Symbol | {__Symbol}] :=
  ($cwNotes = Append[$cwNotes, h]; cw[e]);
cw[D[f_, x_Symbol]] /; valuedQ[x] := cvalue[x, symbolProblem[x]];
cw[D[f_, x_Symbol]] := Module[{fv = cwValue[cw[f]], d},
  If[!FreeQ[fv, cfail | cvalue] || !polyQ[fv, Variables[fv]] || !FreeQ[fv, _Real | _Complex],
    cfail[D],
    d = D[fv, x];
    $cwTrusted = Append[$cwTrusted, "Mathilda's differentiation: D[" <>
      ToString[fv, InputForm] <> ", " <> ToString[x, InputForm] <> "] = " <>
      ToString[d, InputForm] <> " (Lean checks the claim with this derivative written in)"];
    hval[d]]];
cw[(h_Symbol)[a___]] /; MemberQ[$cwKeep, h] := inert[h, cwl[a]];
(* Assumptions may also use || and ! (the identity generator states them). *)
cw[(h:Or | Not)[a___]] /; $cwAssum := inert[h, cwl[a]];
cw[(h_)[___]] := cfail[h];
cw[e_] := cfail[e];
cwl[] := Sequence[];
cwl[a_, r___] := Sequence[cw[a], cwl[r]];

(* A whole side of a claimed equation may be an indefinite Integrate (the
   other side is the claimed antiderivative) or a Sum (the other side is the
   claimed closed form).  Integrate is put on the left. *)
cwClaim[Equal[Integrate[f_, x_Symbol], F_]] :=
  inert[Equal, inert[Integrate, cw[f], csym[x]], Block[{$cwAnti = True}, cw[F]]];
cwClaim[Equal[F_, Integrate[f_, x_Symbol]]] :=
  inert[Equal, inert[Integrate, cw[f], csym[x]], Block[{$cwAnti = True}, cw[F]]];
cwClaim[Equal[Sum[t_, {k_Symbol, a_, n_}], r_]] := Block[{$cwAnyPow = True},
  inert[Equal, inert[Sum, cw[t], inert[List, csym[k], cw[a], cw[n]]], cw[r]]];
cwClaim[(h_Symbol)[a___]] /; MemberQ[$cwComputations, h] && h =!= Sum := inert[h, cwl[a]];
cwClaim[Sum[t_, {k_Symbol, a_, n_}]] := Block[{$cwAnyPow = True},
  inert[Sum, cw[t], inert[List, csym[k], cw[a], cw[n]]]];
cwClaim[e_] := cw[e];

cwHeld[t_] := HoldComplete @@ {t} //. {inert[h_, a___] :> h[a], hk[x_] :> x, hval[v_] :> v};
cwValue[t_] := If[FreeQ[t, cfail | cvalue], ReleaseHold[cwHeld[t]], t];

(* {HoldComplete[rewritten input], stripped heads, trusted steps}, or
   $Failed after Certify::eval / Certify::value.  A symbol holding the input
   (c = a == b; Certify[c]) is refused like any other symbol with a value:
   its value was computed when it was assigned, so it is not the claim as
   written (Certify[Evaluate[c]] certifies the value explicitly). *)
rewriteInput[HoldComplete[c_Symbol]] /; valuedQ[c] :=
  (symbolMessage[HoldForm[c], symbolProblem[c]]; $Failed);
rewriteInput[hc:HoldComplete[c_]] := Module[{t, bad},
  Block[{$cwNotes = {}, $cwTrusted = {}},
    t = cwClaim[c];
    bad = Cases[{t}, cvalue[s_, v_] :> {HoldForm[s], v}, Infinity];
    If[bad =!= {}, symbolMessage @@ First[bad]; Return[$Failed, Module]];
    bad = Cases[{t}, cfail[h_] :> h, Infinity];
    If[bad =!= {},
      Message[Certify::eval, HoldForm[c], First[bad]]; Return[$Failed, Module]];
    {cwHeld[t], $cwNotes, $cwTrusted}]];

(* Mathilda's own evaluation of the two sides of an Equal input (for the
   heuristic engine and "Result"), never used for the Lean statement. *)
mathildaSides[HoldComplete[Equal[a_, b_]]] := {a, b};
mathildaSides[HoldComplete[c_Symbol]] := With[{v = c},
  If[HoldComplete[v] =!= HoldComplete[c], mathildaSides[HoldComplete[v]], None]];
mathildaSides[_] := None;

(* Dispatcher: one rule per supported computation / claim shape. *)
generate[HoldComplete[Equal[a_, b_]], o_] := Module[{av = a, bv = b},
  With[{x = av, y = bv},
    identityGen["Identity", HoldForm[a == b], HoldForm[x == y], av, bv, {}, o,
      {rawH[a], rawH[b]}]]];

generate[HoldComplete[(h_Symbol)[e_, rest___]], o_] /;
    MemberQ[Join[$identityHeads, $factorHeads], h] :=
  Module[{ev = e, res},
    res = h[ev, rest];
    identityGen[
      If[MemberQ[$factorHeads, h], "Factorization", "Identity"],
      HoldForm[h[e, rest]], res, ev, res,
      If[MemberQ[$factorHeads, h], {"irreducibility of the factors (not certified)"}, {}],
      o, {rawH[e], rawFraction[res]}]];

(* ---------------------------------------------------------------------- *)
(* Shared helpers for the structured certificates (rows 4-10)              *)
(* ---------------------------------------------------------------------- *)

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
  "Lean.Private.ratCast " ~~ n:(LetterCharacter ~~ DigitCharacter..) :>
    "(" <> n <> " : Rat)"];
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

(* Inert heads keep a structured expression (p' q - p q', (q^2 + p^2) ...)
   from being normalised; inertHeld turns it into the held statement. *)
$inertRules = {iPlus -> Plus, iTimes -> Times, iPower -> Power, iEqual -> Equal,
  iUnequal -> Unequal, iAnd -> And, iOr -> Or, iImplies -> Implies,
  iEquivalent -> Equivalent, iGreaterEqual -> GreaterEqual, iLessEqual -> LessEqual,
  iNot -> Not, iForAll -> ForAll, iPrimeQ -> PrimeQ, iSumOp -> Sum, iLess -> Less,
  iGreater -> Greater};
inertValue[e_] := e /. $inertRules;
inertHeld[e_] := With[{x = e}, HoldForm[x]] /. $inertRules;

(* iTimes / iPlus that drop unit factors and zero terms. *)
iProd[fs___] := With[{l = DeleteCases[{fs}, 1]},
  Which[l === {}, 1, Length[l] === 1, First[l], True, iTimes @@ l]];
iSum[ts___] := With[{l = DeleteCases[{ts}, 0]},
  Which[l === {}, 0, Length[l] === 1, First[l], True, iPlus @@ l]];

(* ForAll[vars, s] around a held statement s (vars may be empty). *)
forAllHeld[held_HoldForm, {}] := held;
forAllHeld[held_HoldForm, vars_List] := With[{vs = vars},
  held /. HoldForm[s_] :> HoldForm[ForAll[vs, s]]];

ratNumberQ[e_] := MatchQ[e, _Integer | _Rational];
ratMatrixQ[m_] := MatrixQ[m] && AllTrue[Flatten[m], ratNumberQ];
ratVectorQ[v_] := VectorQ[v] && AllTrue[v, ratNumberQ];

(* The generator association every structured kind returns. *)
structGen[kind_, input_, result_, held_, stmt_, proof_, cert_, trusted_, hyps_, claim_,
    dom_:Rationals] := <|"Kind" -> kind, "Input" -> input, "Result" -> result,
  "Theorem" -> held, "Statement" -> stmt, "Proof" -> proof, "Certificate" -> cert,
  "TrustedSteps" -> trusted, "Hypotheses" -> hyps, "Tier" -> 0, "Domain" -> dom,
  "Claim" -> claim|>;

nocert[input_] := (Message[Certify::nocert, input]; $Failed);

(* ---------------------------------------------------------------------- *)
(* Row 4: Integrate[f, x], f in Q(x) -- the Tier 0 cleared chain-rule       *)
(* identity (plan section 3.2)                                             *)
(* ---------------------------------------------------------------------- *)

(* Test-only hook: when not None, this expression replaces the antiderivative
   Integrate computed (it must be in the integration variable).  Used to
   check that Lean rejects a wrong antiderivative. *)
$antiderivOverride = None;
$surdSerial = 0;

(* d^(n/2) with n odd, as d^((n-1)/2) * Sqrt[d], where Sqrt[d] is written over
   the primes: Prod p^Floor[k/2] * Prod_{k odd} sq[p], sq[p] standing for
   Sqrt[p].  Prime generators matter: Q[s_p]/(s_p^2 - p) is the field
   Q(Sqrt[p], ...), so reduction modulo the relations decides zero exactly
   (independent symbols for Sqrt[2], Sqrt[3] and Sqrt[6] would lose
   Sqrt[2] Sqrt[3] = Sqrt[6]). *)
surdPow[1, _] := 1;
surdPow[d_Integer, e_Rational] := d^((Numerator[e] - 1)/2) *
  Times @@ (#[[1]]^Quotient[#[[2]], 2] * If[OddQ[#[[2]]], sq[#[[1]]], 1] & /@
    FactorInteger[d]);

$surdRules = {
  Power[r_Rational?Positive, e_Rational] /; Denominator[e] === 2 :>
    surdPow[Numerator[r], e] * surdPow[Denominator[r], -e],
  Power[d_Integer?Positive, e_Rational] /; Denominator[e] === 2 :> surdPow[d, e]};

(* F is built from the variables by Plus, Times, integer Power and one level
   of Log / ArcTan / ArcTanh. *)
elemQ[e_] := MatchQ[e, (Log | ArcTan | ArcTanh)[_]];
allowedQ[e_, vars_] := Which[
  ratNumberQ[e], True,
  Head[e] === Symbol, MemberQ[vars, e],
  MatchQ[e, _Plus | _Times], AllTrue[List @@ e, allowedQ[#, vars] &],
  MatchQ[e, Power[_, _Integer]], allowedQ[e[[1]], vars],
  elemQ[e], allowedQ[e[[1]], vars] && FreeQ[e[[1]], Log | ArcTan | ArcTanh],
  True, False];

(* The derivative of c * phi(p/q) as {c, N, D} (N/D = d/dx phi(p/q)), with N
   and D in structured, unexpanded form. *)
elemPiece[c_, h_, u_, x_, vars_] := Module[{tu, p, q, pp, qq},
  tu = Together[u]; p = Numerator[tu]; q = Denominator[tu];
  If[!polyQ[p, vars] || !polyQ[q, vars], Return[$Failed, Module]];
  pp = D[p, x]; qq = D[q, x];
  {c, iPlus[iTimes[pp, q], iTimes[-1, p, qq]],
   Switch[h,
     Log, iTimes[p, q],
     ArcTan, iPlus[iPower[q, 2], iPower[p, 2]],
     ArcTanh, iPlus[iPower[q, 2], iTimes[-1, iPower[p, 2]]]],
   h, p, q}];

antiderivGen[input_, f_, x_Symbol, o_, claimed_:None] := Module[{F, tf, a, b, Fs, bases,
    svars, vars, sub, terms, rpart = 0, elems = {}, pieces, tr, r, w, dens, lhs, lhsV, rels,
    red, ks, rem, rhs, held, stmt, nonzero, hyps, back, fs, ph, cf},
  (* The claimed antiderivative of Integrate[f, x] == F, else the test hook,
     else Mathilda's own Integrate. *)
  F = Which[claimed =!= None, claimed, $antiderivOverride =!= None, $antiderivOverride,
    True, Integrate[f, x]];
  If[!FreeQ[F, Integrate], Return[nocert[input], Module]];
  (* The integrand must be a rational function of x over Q. *)
  tf = Together[f]; a = Numerator[tf]; b = Denominator[tf];
  If[!polyQ[a, {x}] || !polyQ[b, {x}] || !SubsetQ[{x}, Variables[{a, b}]] ||
     !FreeQ[{a, b}, _Complex | _Real],
    Return[nocert[input], Module]];
  (* Tier -> 1: the analytic HasDerivAt statement (row 17). *)
  If[o["Tier"] === 1, Return[tier1Or[input, F, f, x, o], Module]];
  (* Quadratic surds Sqrt[d] become fresh variables s_j (relation s_j^2 = d_j). *)
  Fs = F /. $surdRules;
  bases = Union[Cases[Fs, sq[d_] :> d, {0, Infinity}]];
  (* Fresh symbols per call (sixth review): Lean`Private`s<j>$<serial>,
     printed as s<j> by leanRename.  A claim cannot name them (symbols of
     the Lean`Private` context are refused in claims, Certify::context), and
     being fresh they cannot coincide with anything that exists already. *)
  $surdSerial++;
  svars = Table[Symbol["Lean`Private`" <>
      If[StringMatchQ[SymbolName[x], "s" ~~ DigitCharacter..], "r", "s"] <> ToString[j] <>
      "$" <> ToString[$surdSerial]],
    {j, Length[bases]}];
  sub = Thread[(sq /@ bases) -> svars];
  Fs = Fs /. sub;
  vars = Join[{x}, svars];
  (* Outside the Tier 0 shape: Tier 1 when its HasDerivAt builder applies
     (prepare then needs a project), Certify::tier otherwise. *)
  If[!allowedQ[Fs, vars], Return[tier1Or[input, F, f, x, o], Module]];
  (* Decompose F = R + Sum c_i phi_i(u_i), each c_i free of x. *)
  terms = If[Head[Fs] === Plus, List @@ Fs, {Fs}];
  Do[
    If[FreeQ[t, Log | ArcTan | ArcTanh], rpart += t,
      fs = If[Head[t] === Times, List @@ t, {t}];
      ph = Select[fs, elemQ]; cf = Times @@ Select[fs, !elemQ[#] &];
      If[Length[ph] =!= 1 || !FreeQ[cf, x], Return[nocert[input], Module]];
      AppendTo[elems, {cf, Head[First[ph]], First[First[ph]]}]],
    {t, terms}];
  pieces = elemPiece[#[[1]], #[[2]], #[[3]], x, vars] & /@ elems;
  If[MemberQ[pieces, $Failed], Return[nocert[input], Module]];
  tr = Together[rpart]; r = Numerator[tr]; w = Denominator[tr];
  If[!polyQ[r, vars] || !polyQ[w, vars], Return[nocert[input], Module]];
  If[r =!= 0, AppendTo[pieces, If[FreeQ[w, x],
    {1, D[r, x], w, "Rational", r, w},
    {1, iPlus[iTimes[D[r, x], w], iTimes[-1, r, D[w, x]]], iPower[w, 2], "Rational", r, w}]]];
  (* Structured cleared identity:
       (Sum_i c_i N_i Prod_{j != i} D_j) b - a Prod_j D_j = Sum_j (s_j^2 - d_j) K_j *)
  dens = pieces[[All, 3]];
  lhs = iSum[Sequence @@ Table[
      iProd[pieces[[i, 1]], pieces[[i, 2]], Sequence @@ Delete[dens, i], b],
      {i, Length[pieces]}],
    iProd[-1, a, Sequence @@ dens]];
  rels = MapThread[#1^2 - #2 &, {svars, bases}];
  red[e_] := If[svars === {}, {{}, Expand[e]},
    PolynomialReduce[Expand[e], rels, vars]];
  lhsV = Expand[inertValue[lhs]];
  {ks, rem} = red[lhsV];
  (* Mathilda-side checks: a denominator that vanishes identically would make
     the cleared identity vacuous, so there is no certificate. *)
  back = Thread[svars -> (Sqrt /@ bases)];
  nonzero[e_] := Last[red[inertValue[e]]] =!= 0 &&
    Expand[inertValue[e] /. back] =!= 0;
  If[!AllTrue[dens, nonzero], Return[nocert[input], Module]];
  rhs = iSum @@ MapThread[If[#3 === 0, 0, iTimes[iPlus[iPower[#1, 2], -#2], #3]] &,
    {svars, bases, ks}];
  held = forAllHeld[inertHeld[iEqual[lhs, rhs]], vars];
  stmt = leanThm[held, "Rat"];
  If[!StringQ[stmt], Return[nocert[input], Module]];
  hyps = DeleteDuplicates[Cases[Append[dens, b],
    d_ /; !NumberQ[inertValue[d]] :> (inertHeld[iUnequal[d, 0]] /. back)]];
  structGen["Antiderivative", input, F, held, stmt,
    "by" <> introLine[vars] <> "\n  grind",
    <|"Form" -> "ClearedDerivative", "Variable" -> x, "Integrand" -> {a, b},
      "Terms" -> Map[<|"Coefficient" -> #[[1]], "Function" -> #[[4]],
          "Numerator" -> #[[5]], "Denominator" -> #[[6]]|> &, pieces],
      "Surds" -> Thread[svars -> (Sqrt /@ bases)], "SurdCofactors" -> ks,
      "PreCheck" -> (rem === 0)|>,
    {"elementary derivative rules d log u = u'/u, d arctan u = u'/(1 + u^2), \
d artanh u = u'/(1 - u^2) and the quotient rule (applied symbolically by Mathilda; Lean \
checks only the resulting polynomial identity)",
     "polynomial differentiation by Mathilda's D (the derivatives written into the \
statement)",
     "Together normal form of the integrand (f = a/b), of each argument u_i = p_i/q_i \
and of the rational part R = r/w of the antiderivative",
     "the cleared denominators are not identically zero (checked by Mathilda; listed in \
\"Hypotheses\")"},
    hyps, {D[F, x], f}]];

(* The Tier 1 generator, or $Failed after Certify::nocert (Tier -> 1
   requested) / Certify::tier (no certificate at either tier). *)
tier1Or[input_, F_, f_, x_, o_] := With[{g = tier1AntiderivGen[input, F, f, x, o]},
  If[AssociationQ[g], g,
    If[o["Tier"] === 1, Message[Certify::nocert, input], Message[Certify::tier, input]];
    $Failed]];

generate[HoldComplete[Integrate[f_, x_Symbol]], o_] :=
  antiderivGen[HoldForm[Integrate[f, x]], f, x, o];

(* The claim Integrate[f, x] == F (rewriteInput puts Integrate on the left):
   F, the user's antiderivative, is checked -- not Mathilda's.  The claim as
   written is NOT what is certified (sixth review, v0.6522): in Mathematica
   Integrate[2 x, x] == x^2 + 5 is False, since Integrate returns one
   particular antiderivative, and Lean cannot decide which one without
   evaluating Integrate.  What Lean checks is that F is an antiderivative of
   f, i.e. the claim D[F, x] == f; that is what "Theorem" shows, with
   "Input" still the claim as written, a Certify::antideriv message, and a
   trusted step that states the gap: Integrate[f, x] == F holds only up to
   an additive term whose derivative is 0 -- locally constant, i.e.
   constant on each interval where F and f are defined, and possibly
   different across a singularity (ArcTan[1/x] and -ArcTan[x] differ by
   Pi/2 for x > 0 and by -Pi/2 for x < 0). *)
generate[HoldComplete[Equal[Integrate[f_, x_Symbol], F_]], o_] := Module[{g},
  g = antiderivGen[HoldForm[Integrate[f, x] == F], f, x, o, F];
  If[!AssociationQ[g], Return[g, Module]];
  Message[Certify::antideriv, HoldForm[Integrate[f, x] == F], HoldForm[D[F, x] == f]];
  Join[g, <|
    "Theorem" -> If[g["Tier"] === 1, g["Theorem"], HoldForm[ForAll[x, D[F, x] == f]]],
    "TrustedSteps" -> Append[g["TrustedSteps"],
      "the claim Integrate[f, x] == F is certified as D[F, x] == f (F is an antiderivative \
of f, shown in \"Theorem\"); Integrate[f, x] == F itself holds only up to an additive term \
with derivative 0, which is locally constant (constant on each interval where F and f are \
defined; it may differ across a point where they are not, as ArcTan[1/x] and -ArcTan[x] \
do at x = 0), not checked"]|>]];

(* ---------------------------------------------------------------------- *)
(* Row 5: PolynomialExtendedGCD (Bezout)                                   *)
(* ---------------------------------------------------------------------- *)

generate[HoldComplete[PolynomialExtendedGCD[a_, b_, x_Symbol]], o_] :=
  Module[{input = HoldForm[PolynomialExtendedGCD[a, b, x]], av = a, bv = b, res, g, s, t,
      qa, qb, held, ty, stmt},
  If[!polyQ[av, {x}] || !polyQ[bv, {x}] || !SubsetQ[{x}, Variables[{av, bv}]],
    Return[nocert[input], Module]];
  res = PolynomialExtendedGCD[av, bv, x];
  If[!MatchQ[res, {_, {_, _}}], Return[nocert[input], Module]];
  {g, {s, t}} = res;
  qa = PolynomialQuotient[av, g, x]; qb = PolynomialQuotient[bv, g, x];
  held = forAllHeld[inertHeld[iAnd[
      iEqual[iPlus[iTimes[s, av], iTimes[t, bv]], g],
      iEqual[av, iTimes[g, qa]], iEqual[bv, iTimes[g, qb]]]], {x}];
  ty = autoDomain[{av, bv, g, s, t, qa, qb}];
  stmt = leanThm[held, ty];
  If[!StringQ[stmt], Return[nocert[input], Module]];
  structGen["Bezout", input, res, held, stmt,
    "by" <> introLine[{x}] <> "\n  refine And.intro ?_ (And.intro ?_ ?_) <;> grind",
    <|"GCD" -> g, "Cofactors" -> {s, t}, "Quotients" -> {qa, qb},
      "PreCheck" -> (Expand[s av + t bv - g] === 0 && Expand[av - g qa] === 0 &&
        Expand[bv - g qb] === 0)|>,
    {}, {}, {{s av + t bv, av, bv}, {g, g qa, g qb}}, leanType[ty]]];

(* ---------------------------------------------------------------------- *)
(* Row 6: PolynomialReduce and ideal membership Implies[p1 == 0 && ..., q == 0] *)
(* ---------------------------------------------------------------------- *)

monomials[vars_List, d_Integer] := Times @@ (vars^#) & /@
  Select[Tuples[Range[0, d], Length[vars]], Total[#] <= d &];

(* Cofactors c with q == c . ps: PolynomialReduce first, then a degree-bounded
   ansatz c_i = Sum_{|alpha| <= d} u_{i,alpha} x^alpha (d = 0..dmax) whose
   coefficient equations are solved with LinearSolve (free parameters 0). *)
cofactors[q_, ps_List, vars_List, dmax_Integer:4] := Module[{pr, monos, cols, rules,
    qr, keys, vec, sol, nm, cs},
  pr = Quiet[PolynomialReduce[q, ps, vars]];
  If[MatchQ[pr, {_List, _}] && Expand[q - pr[[1]] . ps] === 0, Return[pr[[1]], Module]];
  vec[rl_] := Replace[#, Append[rl, _ -> 0]] & /@ keys;
  Do[
    monos = monomials[vars, d]; nm = Length[monos];
    cols = Flatten[Table[Expand[m p], {p, ps}, {m, monos}], 1];
    rules = CoefficientRules[#, vars] & /@ cols;
    qr = CoefficientRules[Expand[q], vars];
    keys = Union @@ Join[(First /@ #) & /@ rules, {First /@ qr}];
    sol = Quiet[LinearSolve[Transpose[vec /@ rules], vec[qr]]];
    If[VectorQ[sol],
      cs = Table[sol[[(i - 1) nm + 1 ;; i nm]] . monos, {i, Length[ps]}];
      If[Expand[q - cs . ps] === 0, Return[cs, Module]]],
    {d, 0, dmax}];
  $Failed];

generate[HoldComplete[PolynomialReduce[q_, ps_List, vs_List]], o_] :=
  Module[{input = HoldForm[PolynomialReduce[q, ps, vs]], qv = q, psv = ps, res, cs, r,
      vars, held, ty, stmt},
  vars = Sort[Union[Variables[{qv, psv}], vs]];
  If[!AllTrue[vars, Head[#] === Symbol &] || !polyQ[qv, vars] ||
     !AllTrue[psv, polyQ[#, vars] &], Return[nocert[input], Module]];
  res = PolynomialReduce[qv, psv, vs];
  If[!MatchQ[res, {_List, _}], Return[nocert[input], Module]];
  {cs, r} = res;
  held = forAllHeld[inertHeld[iEqual[qv,
      iSum[Sequence @@ MapThread[iTimes, {cs, psv}], r]]], vars];
  ty = autoDomain[{qv, psv, cs, r}];
  stmt = leanThm[held, ty];
  If[!StringQ[stmt], Return[nocert[input], Module]];
  structGen["IdealMembership", input, res, held, stmt,
    "by" <> introLine[vars] <> "\n  grind",
    <|"Cofactors" -> cs, "Remainder" -> r, "Divisors" -> psv,
      "PreCheck" -> (Expand[qv - cs . psv - r] === 0)|>,
    {"that the remainder is fully reduced with respect to the divisors (not certified; \
the division identity is)"},
    {}, {qv, cs . psv + r}, leanType[ty]]];

(* Hypotheses l_i == r_i, goal l == r: q = l - r lies in the ideal of the
   p_i = l_i - r_i.  Proof shape (plan section 0.1): have k_i : p_i = 0,
   have hc : q = Sum c_i * p_i, rw [k_i] at hc, grind. *)
membershipGen[input_, hs_List, {l_, r_}, o_] := Module[{ps, q, vars, cs, eq, held, ty,
    stmt, ptxt, qtxt, ctxt, proof, k},
  ps = Expand[#[[1]] - #[[2]]] & /@ hs; q = Expand[l - r];
  vars = Sort[Variables[{ps, q}]];
  If[!AllTrue[vars, Head[#] === Symbol &] || !AllTrue[Append[ps, q], polyQ[#, vars] &],
    Return[nocert[input], Module]];
  cs = cofactors[q, ps, vars];
  If[cs === $Failed, Return[nocert[input], Module]];
  eq[{a_, b_}] := iEqual[a, b];
  held = forAllHeld[inertHeld[Fold[iImplies[#2, #1] &, eq[{l, r}], Reverse[eq /@ hs]]],
    vars];
  (* Hypotheses: over Int (or Q) the implication would be weaker than
     intended (fewer points satisfy them), so the statement is over every
     field of characteristic 0; the cofactor identity proved inside it holds
     in every commutative ring. *)
  ty = "Rat";
  stmt = leanThm[held, ty];
  ptxt = leanStmt[inertHeld[iPlus[#[[1]], iTimes[-1, #[[2]]]]], ty] & /@ hs;
  qtxt = leanStmt[inertHeld[iPlus[l, iTimes[-1, r]]], ty];
  ctxt = leanStmt[HoldForm[#], ty] & /@ cs;
  If[!StringQ[stmt] || !AllTrue[Join[ptxt, ctxt, {qtxt}], StringQ],
    Return[nocert[input], Module]];
  k = Length[ps];
  proof = If[MemberQ[ptxt, qtxt], "by" <> introLine[vars] <> "\n  intros\n  grind",
    StringJoin["by", introLine[vars], "\n  intros",
      Table["\n  have mth_k" <> ToString[i] <> " : " <> ptxt[[i]] <> " = 0 := (by grind)",
        {i, k}],
      "\n  have mth_hc : ", qtxt, " = ",
      StringRiffle[Table["(" <> ctxt[[i]] <> ") * (" <> ptxt[[i]] <> ")", {i, k}], " + "],
      " := (by grind)",
      "\n  rw [", StringRiffle[Table["mth_k" <> ToString[i], {i, k}], ", "], "] at mth_hc",
      "\n  grind"]];
  structGen["IdealMembership", input, cs, held, genericStmt[stmt, "CharZero"],
    genericText[proof], <|"Cofactors" -> cs, "Generators" -> ps, "Goal" -> q|>, {}, {},
    {q, cs . ps}, $genericDomain["CharZero"]]];

generate[HoldComplete[Implies[And[hs__Equal], Equal[l_, r_]]], o_] :=
  membershipGen[HoldForm[Implies[And[hs], l == r]],
    Cases[HoldComplete[hs], HoldPattern[Equal[a_, b_]] :> {a, b}], {l, r}, o];

generate[HoldComplete[Implies[h_Equal, Equal[l_, r_]]], o_] :=
  membershipGen[HoldForm[Implies[h, l == r]],
    Cases[HoldComplete[h], HoldPattern[Equal[a_, b_]] :> {a, b}], {l, r}, o];

(* ---------------------------------------------------------------------- *)
(* Row 7: Solve -- completeness for univariate polynomials with rational   *)
(* roots, soundness (each solution satisfies each equation) otherwise      *)
(* ---------------------------------------------------------------------- *)

(* Multiplicity of the rational root rt of p. *)
rootMultiplicity[p_, x_, rt_] := Module[{m = 0, f = p},
  While[f =!= 0 && Expand[PolynomialRemainder[f, x - rt, x]] === 0,
    m++; f = PolynomialQuotient[f, x - rt, x]];
  m];

(* The relations Solve accepts in its first argument, with their inert heads. *)
$solveRels = {Equal -> iEqual, Unequal -> iUnequal, Less -> iLess, LessEqual -> iLessEqual,
  Greater -> iGreater, GreaterEqual -> iGreaterEqual};

(* The conjuncts of a held Solve system (And / List nest flattened). *)
solveConjuncts[HoldComplete[And[a___]]] := Flatten[solveConjuncts /@ Thread[HoldComplete[{a}]], 1];
solveConjuncts[HoldComplete[List[a___]]] := Flatten[solveConjuncts /@ Thread[HoldComplete[{a}]], 1];
solveConjuncts[h_HoldComplete] := {h};

(* One conjunct as {head, lhs, rhs} with evaluated sides; $Failed when it is
   not a binary relation of $solveRels (Or, Element, chains: no certificate). *)
solveConstraint[HoldComplete[(h_Symbol)[a_, b_]]] /; MemberQ[Keys[$solveRels], h] :=
  {h, a, b};
solveConstraint[_] := $Failed;

(* Bases of the negative powers (denominators) in e. *)
denominatorBases[e_] := DeleteDuplicates[Cases[{e}, Power[b_, _?Negative] :> b,
  {0, Infinity}]];

solveGen[input_, cons_List, vars_List, sol_, o_] := Module[{x, p, roots, lc, mult, fact,
    held, stmt, ftxt, etxt, l, r, alts, sound, dens, eqs},
  If[!MatchQ[sol, {{__Rule} ..}] || !AllTrue[Flatten[vars /. sol], ratNumberQ],
    Return[nocert[input], Module]];
  eqs = Select[cons, #[[1]] === Equal &][[All, 2 ;; 3]];
  (* Completeness: the system is one univariate polynomial equation (no
     other constraint) that splits over Q. *)
  If[Length[vars] === 1 && Length[cons] === 1 && Length[eqs] === 1,
    x = First[vars]; {l, r} = First[eqs]; p = Expand[l - r];
    If[polyQ[p, {x}] && polyQ[l, {x}] && polyQ[r, {x}] && SubsetQ[{x}, Variables[p]] &&
       Exponent[p, x] > 0,
      roots = DeleteDuplicates[x /. sol];
      lc = Coefficient[p, x, Exponent[p, x]];
      mult = rootMultiplicity[p, x, #] & /@ roots;
      fact = iProd[lc, Sequence @@ MapThread[
        If[#2 === 1, iPlus[x, -#1], iPower[iPlus[x, -#1], #2]] &, {roots, mult}]];
      If[Expand[p - inertValue[fact]] === 0,
        alts = iEqual[x, #] & /@ roots;
        held = forAllHeld[inertHeld[iEquivalent[iEqual[l, r],
            If[Length[alts] === 1, First[alts], iOr @@ alts]]], {x}];
        stmt = leanThm[held, "Rat"];
        ftxt = leanStmt[inertHeld[fact], "Rat"];
        etxt = leanStmt[inertHeld[iEqual[l, r]], "Rat"];
        If[StringQ[stmt] && StringQ[ftxt] && StringQ[etxt],
          Return[structGen["Roots", input, sol, held, stmt,
            StringJoin["by", introLine[{x}],
              "\n  have mth_e : ", leanStmt[inertHeld[iPlus[l, iTimes[-1, r]]], "Rat"],
              " = ", ftxt, " := (by grind)",
              "\n  have mth_h : (", etxt, ") \[LeftRightArrow] ", ftxt,
              " = 0 := (by rw [\[LeftArrow] mth_e]; grind)",
              "\n  rw [mth_h]\n  grind"],
            <|"Completeness" -> True, "Factorization" -> inertValue[fact],
              "Roots" -> roots, "Multiplicities" -> mult|>,
            {}, {}, {Flatten[(l - r) /. sol], 0}], Module]]]]];
  (* Soundness: substitute every solution into every constraint (held, so
     the substitution reaches Lean) and let the kernel evaluate.  Lean's
     x / 0 = 0 would accept a "solution" that zeroes a denominator, so each
     denominator of the system is also stated nonzero at each solution. *)
  dens = denominatorBases[cons[[All, 2 ;; 3]]];
  sound = With[{e = iAnd @@ Join[iUnequal[#, 0] & /@ dens,
       (#[[1]] /. $solveRels) @@ #[[2 ;; 3]] & /@ cons]},
    HoldForm[e]] /. iAnd[c_] :> c;
  held = inertHeld[If[Length[sol] === 1, sound /. First[sol], iAnd @@ ((sound /. #) & /@ sol)]];
  stmt = leanThm[held, "Rat"];
  If[!StringQ[stmt], Return[nocert[input], Module]];
  structGen["Roots", input, sol, held, stmt, "by\n  decide +kernel",
    <|"Completeness" -> False, "Solutions" -> sol, "Denominators" -> dens|>,
    {"completeness: that Solve returned every solution (only soundness -- each solution \
satisfies each constraint, with every denominator nonzero -- is certified)"},
    {}, {Flatten[Table[(#[[1]] - #[[2]]) & /@ eqs /. s, {s, sol}]], 0}]];

generate[HoldComplete[Solve[eqs_, vs_]], o_] := Module[{input = HoldForm[Solve[eqs, vs]],
    cons, vars, sol},
  cons = solveConstraint /@ solveConjuncts[HoldComplete[eqs]];
  vars = Flatten[{vs}];
  If[cons === {} || MemberQ[cons, $Failed] || !MemberQ[cons[[All, 1]], Equal] ||
     !AllTrue[vars, Head[#] === Symbol &], Return[nocert[input], Module]];
  sol = Quiet[Solve[eqs, vs]];
  solveGen[input, cons, vars, sol, o]];

(* ---------------------------------------------------------------------- *)
(* Rows 8-10: Inverse, Det (n <= 6), LinearSolve, NullSpace over Q,        *)
(* checked by the kernel evaluating the prelude's matMul / det / matVec    *)
(* ---------------------------------------------------------------------- *)

matrixGen[kind_, input_, result_, held_, cert_, trusted_, claim_] := With[
  {stmt = leanThm[held, "Rat"]},
  If[!StringQ[stmt], nocert[input],
    structGen[kind, input, result, held, stmt, "by\n  decide +kernel", cert, trusted, {},
      claim]]];

generate[HoldComplete[Inverse[m_]], o_] := Module[{input = HoldForm[Inverse[m]], a = m,
    b, n},
  If[!ratMatrixQ[a] || Length[a] =!= Length[First[a]], Return[nocert[input], Module]];
  b = Quiet[Inverse[a]];
  If[!ratMatrixQ[b], Return[nocert[input], Module]];
  n = Length[a];
  If[Dimensions[b] =!= {n, n}, Return[nocert[input], Module]];
  matrixGen["MatrixInverse", input, b,
    With[{aa = a, bb = b, k = n},
      HoldForm[And[hasShape[aa, k, k], hasShape[bb, k, k],
        matMul[aa, bb] == idMat[k], matMul[bb, aa] == idMat[k]]]],
    <|"Matrix" -> a, "Inverse" -> b|>, {}, {{a . b, b . a}, {IdentityMatrix[n], IdentityMatrix[n]}}]];

generate[HoldComplete[Det[m_]], o_] := Module[{input = HoldForm[Det[m]], a = m, d},
  If[!ratMatrixQ[a] || Length[a] =!= Length[First[a]] || Length[a] > 6,
    Return[nocert[input], Module]];
  d = Det[a];
  matrixGen["Determinant", input, d, With[{aa = a, dd = d, k = Length[a]},
      HoldForm[hasShape[aa, k, k] && det[aa] == dd]],
    <|"Matrix" -> a, "Method" -> "Laplace expansion along row 0 (prelude Mathilda.det)"|>,
    {}, {d, Det[a]}]];

generate[HoldComplete[LinearSolve[m_, v_]], o_] := Module[{input = HoldForm[LinearSolve[m, v]],
    a = m, b = v, x},
  If[!ratMatrixQ[a] || !ratVectorQ[b] || Length[b] =!= Length[a],
    Return[nocert[input], Module]];
  x = Quiet[LinearSolve[a, b]];
  If[!ratVectorQ[x], Return[nocert[input], Module]];
  matrixGen["LinearSolve", input, x, With[{aa = a, xx = x, bb = b, r = Length[a],
      c = Length[First[a]]},
      HoldForm[And[hasShape[aa, r, c], hasLength[xx, c], matVec[aa, xx] == bb]]],
    <|"Matrix" -> a, "RightHandSide" -> b, "Solution" -> x|>, {}, {a . x, b}]];

generate[HoldComplete[NullSpace[m_]], o_] := Module[{input = HoldForm[NullSpace[m]], a = m,
    vs, zero},
  If[!ratMatrixQ[a], Return[nocert[input], Module]];
  vs = NullSpace[a];
  If[!MatchQ[vs, {__}] || !AllTrue[vs, ratVectorQ], Return[nocert[input], Module]];
  zero = ConstantArray[0, Length[a]];
  matrixGen["NullSpace", input, vs,
    inertHeld[iAnd @@ Join[{hasShape[a, Length[a], Length[First[a]]]},
      Flatten[With[{aa = a, z = zero, c = Length[First[a]]},
        {hasLength[#, c], iEqual[matVec[aa, #], z]}] & /@ vs, 1]]],
    <|"Matrix" -> a, "Basis" -> vs|>,
    {"completeness: that the vectors are independent and span the whole null space \
(only A.v = 0 for each v is certified)"},
    {a . Transpose[vs], ConstantArray[0, {Length[a], Length[vs]}]}]];

(* ---------------------------------------------------------------------- *)
(* Rows 12-13: sums of squares -- ForAll[vars, p >= q] with optional        *)
(* g_i >= 0 hypotheses (Putinar).  Lean gets the hints 0 <= l_k^2 and       *)
(* 0 <= l_k^2 g, the explicit identity p - q = Sum d_k l_k^2 (g), and grind *)
(* ---------------------------------------------------------------------- *)

(* Total degree of the polynomial p in vars. *)
totalDegree[p_, vars_] := With[{cr = CoefficientRules[Expand[p], vars]},
  If[cr === {} || cr === 0, 0, Max[Total /@ cr[[All, 1]]]]];

(* Coefficient of the monomial m (in vars) in the expanded polynomial p. *)
monoCoeff[rules_, m_, vars_] := Replace[Exponent[m, #] & /@ vars, Append[rules, _ -> 0]];

(* Exact rational LDL^T of a symmetric matrix: {L, d} with q == L.DiagonalMatrix[d].L^T,
   L unit lower triangular; $Failed when q is not positive semidefinite.  No
   pivoting is needed: in a PSD matrix (and in each of its Schur complements)
   a zero diagonal entry has a zero row, which is exactly what is checked. *)
ldl[q_] := Module[{n, L, d, v},
  If[!ratMatrixQ[q] || Length[q] =!= Length[First[q]] || q =!= Transpose[q],
    Return[$Failed, Module]];
  n = Length[q]; L = IdentityMatrix[n]; d = ConstantArray[0, n];
  Do[
    d[[k]] = q[[k, k]] - Sum[L[[k, j]]^2 d[[j]], {j, 1, k - 1}];
    If[d[[k]] < 0, Return[$Failed, Module]];
    Do[
      v = q[[i, k]] - Sum[L[[i, j]] L[[k, j]] d[[j]], {j, 1, k - 1}];
      If[d[[k]] === 0, If[v =!= 0, Return[$Failed, Module]], L[[i, k]] = v/d[[k]]],
      {i, k + 1, n}],
    {k, n}];
  {L, d}];

(* The squares of one Gram block: {{d_k, l_k, mult}, ...} with
   basis.Q.basis == Sum d_k l_k^2 and l = L^T.basis; $Failed if Q is not PSD. *)
gramTerms[mult_, basis_List, q_] := With[{r = ldl[q]},
  If[r === $Failed, $Failed,
    Select[Table[{r[[2, k]], Expand[Sum[r[[1, i, k]] basis[[i]], {i, Length[basis]}]], mult},
      {k, Length[basis]}], #[[1]] =!= 0 &]]];

(* A quadratic form's Gram matrix over {1, x1, ..., xn} is unique. *)
quadGram[p_, vars_] := Module[{z = Prepend[vars, 1], rules = CoefficientRules[Expand[p], vars]},
  Table[If[i === j, monoCoeff[rules, z[[i]]^2, vars], monoCoeff[rules, z[[i]] z[[j]], vars]/2],
    {i, Length[z]}, {j, Length[z]}]];

(* Exact Gram search for q >= 0 (no hypotheses) when the numeric engine
   declines: the Gram matrix over the half-degree monomials is pinned by
   coefficient matching, by Q.m(z) = 0 at every rational real zero z of q
   (a PSD Q with z^T Q z = 0 annihilates m(z)), and by zero rows wherever a
   diagonal entry is forced to 0.  Free parameters are then set to 0 and
   LDL^T decides.  Handles the not-strictly-feasible SDPs (several zeros,
   e.g. (x^2 - 1)^2 + y^2) that a numeric interior point cannot round. *)
exactGramTerms[p_, vars_] := Module[{deg, maxe, basis, n, cs, qs, eqs, rules, crit, zeros,
    sol, qv, zrow, free, prev},
  deg = totalDegree[p, vars];
  If[OddQ[deg] || deg > 8 || Length[vars] > 4, Return[$Failed, Module]];
  maxe = Exponent[p, #] & /@ vars;
  basis = Select[monomials[vars, deg/2],
    Function[m, AllTrue[Range[Length[vars]], Exponent[m, vars[[#]]] <= maxe[[#]]/2 &]]];
  n = Length[basis];
  cs = Table[Unique["Lean`Private`g"], {n (n + 1)/2}];
  qs = Table[With[{a = Min[i, j], b = Max[i, j]}, cs[[(a - 1) n - (a - 1) (a - 2)/2 + (b - a) + 1]]],
    {i, n}, {j, n}];
  rules = CoefficientRules[Expand[basis . qs . basis - p], vars];
  eqs = (#[[2]] == 0) & /@ rules;
  crit = Quiet[Solve[Thread[(D[p, #] & /@ vars) == 0], vars]];
  zeros = If[MatchQ[crit, {{__Rule} ..}],
    Select[crit, AllTrue[vars /. #, ratNumberQ] && Expand[p /. #] === 0 &], {}];
  eqs = Join[eqs, Flatten[Thread[(qs . (basis /. #)) == 0] & /@ zeros]];
  prev = None;
  While[True,
    sol = Quiet[Solve[eqs, cs]];
    If[!MatchQ[sol, {{___Rule}, ___}], Return[$Failed, Module]];
    qv = qs /. First[sol];
    zrow = Select[Range[n], qv[[#, #]] === 0 &];
    If[zrow === prev, Break[]];
    prev = zrow;
    eqs = Join[eqs, Flatten[Thread[qs[[#]] == 0] & /@ zrow]]];
  free = Union[Cases[qv, _Symbol?(MemberQ[cs, #] &), {0, Infinity}]];
  qv = qv /. Thread[free -> 0];
  gramTerms[1, basis, qv]];

(* Test hook: the number of Reduce`SOSCertificate calls (the quadratic path
   never makes one). *)
$sosCalls = 0;
sosCertificate[p_, gs_, vars_] := ($sosCalls++; Reduce`SOSCertificate[p, gs, vars]);

(* {{d, l, mult}, ...} with p == Sum d l^2 mult identically (every d > 0 and
   every mult 1 or one of gs), plus the method used; $Failed if none. *)
sosTerms[p_, gs_List, vars_List] := Module[{q = Expand[p], terms = $Failed, method, cert},
  If[gs === {} && totalDegree[q, vars] <= 2,
    method = "Quadratic"; terms = gramTerms[1, Prepend[vars, 1], quadGram[q, vars]],
    cert = sosCertificate[q, gs, vars];
    If[ListQ[cert],
      method = "SOSCertificate";
      terms = gramTerms @@@ cert;
      terms = If[MemberQ[terms, $Failed], $Failed, Flatten[terms, 1]]];
    If[terms === $Failed && gs === {},
      method = "ExactGram"; terms = exactGramTerms[q, vars]]];
  If[ListQ[terms] && Expand[q - Total[(#[[1]] #[[2]]^2 #[[3]]) & /@ terms]] === 0,
    {terms, method}, $Failed]];

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

(* Hypotheses g >= 0 (or >, <=, <) as polynomials g; $Failed otherwise. *)
ineqPoly[GreaterEqual[a_, b_] | Greater[a_, b_]] := a - b;
ineqPoly[LessEqual[a_, b_] | Less[a_, b_]] := b - a;
ineqPoly[_] := $Failed;
condList[True] := {};
condList[c_And] := List @@ c;
condList[c_] := {c};

(* The relation p >= q (GreaterEqual) or q <= p (LessEqual) held inert. *)
inertRel[GreaterEqual, a_, b_] := iGreaterEqual[a, b];
inertRel[LessEqual, a_, b_] := iLessEqual[a, b];

sosClaimGen[input_, vars_List, conds_List, rel_, a_, b_, o_] := Module[{p, q, gs, target,
    st, terms, gnames, held, stmt, proof, ghints, used},
  If[!AllTrue[vars, Head[#] === Symbol &] || vars === {}, Return[nocert[input], Module]];
  {p, q} = If[rel === GreaterEqual, {a, b}, {b, a}];
  gs = Expand[ineqPoly /@ conds];
  target = Expand[p - q];
  If[MemberQ[gs, $Failed] || !AllTrue[Append[gs, target], polyQ[#, vars] &] ||
     !SubsetQ[vars, Variables[Append[gs, target]]] || !FreeQ[{gs, target}, _Real | _Complex],
    Return[nocert[input], Module]];
  st = sosTerms[target, DeleteCases[gs, _?NumberQ], vars];
  If[st === $Failed, Return[nocert[input], Module]];
  terms = st[[1]];
  used = DeleteCases[Union[terms[[All, 3]]], 1];
  gnames = Table[used[[i]] -> "mth_g" <> ToString[i], {i, Length[used]}];
  held = inertHeld[If[conds === {}, iForAll[vars, inertRel[rel, a, b]],
    iForAll[vars, If[Length[conds] === 1, First[conds], iAnd @@ conds], inertRel[rel, a, b]]]];
  stmt = leanThm[held, "Rat"];
  If[!StringQ[stmt], Return[nocert[input], Module]];
  ghints = Table["\n  have mth_g" <> ToString[i] <> " : 0 \[LessEqual] " <> ratText[used[[i]]] <>
    " := (by grind)", {i, Length[used]}];
  proof = StringJoin["by", introLine[vars], If[conds === {}, "", "\n  intros"], ghints,
    sqHints[terms, gnames],
    "\n  have mth_e : ", ratText[target], " = ", If[terms === {}, "0", sqSum[terms]],
    " := (by grind)", "\n  grind"];
  structGen["Nonnegativity", input, True, held, genericStmt[stmt, "Ordered"],
    genericText[proof],
    <|"Squares" -> terms, "Method" -> st[[2]], "Target" -> target, "Multipliers" -> used,
      "PreCheck" -> True|>,
    {}, {}, None, $genericDomain["Ordered"]]];

generate[HoldComplete[ForAll[vs_, (r:GreaterEqual | LessEqual)[a_, b_]]], o_] :=
  sosClaimGen[HoldForm[ForAll[vs, r[a, b]]], Flatten[{vs}], {}, r, a, b, o];

generate[HoldComplete[ForAll[vs_, c_, (r:GreaterEqual | LessEqual)[a_, b_]]], o_] :=
  sosClaimGen[HoldForm[ForAll[vs, c, r[a, b]]], Flatten[{vs}], condList[c], r, a, b, o];

generate[HoldComplete[Implies[c_, (r:GreaterEqual | LessEqual)[a_, b_]]], o_] := With[{cv = c, av = a, bv = b},
  sosClaimGen[HoldForm[Implies[c, r[a, b]]], Sort[Variables[{ineqPoly /@ condList[cv], av - bv}]],
    condList[cv], r, av, bv, o]];

generate[HoldComplete[(r:GreaterEqual | LessEqual)[a_, b_]], o_] := With[{av = a, bv = b},
  sosClaimGen[HoldForm[r[a, b]], Sort[Variables[av - bv]], {}, r, av, bv, o]];

(* ---------------------------------------------------------------------- *)
(* Row 14: Minimize / Maximize -- KKT multipliers lambda >= 0 at x* and     *)
(* f - m - Sum lambda_j g_j a sum of squares.  The theorem is              *)
(*   (forall x, cons -> f >= m) and cons(xs) and f(xs) = m, xs the optimum *)
(* (<= for Maximize).  Linear equalities are eliminated by subst.  When    *)
(* sigma = 0 (an LP) the multipliers are the dual certificate: weak duality *)
(* plus attainment is optimality.                                          *)
(* ---------------------------------------------------------------------- *)

(* An evaluated polynomial as an inert tree, so that substituting numbers
   keeps the structure (1/2 + 1/2 >= 1, not True). *)
inertOf[e_] := ReleaseHold[With[{v = e}, Hold[v]] /. {Plus -> iPlus, Times -> iTimes,
  Power -> iPower}];

(* Constraints as {kind, poly, inert relation}: "ineq" g >= 0, "eq" h == 0. *)
optCons[c_] := Module[{l = condList[c] /. {LessEqual[x_, y_, z_] :> Sequence[x <= y, y <= z],
    GreaterEqual[x_, y_, z_] :> Sequence[x >= y, y >= z]}},
  (* Heads, not Equal[_, _]: that pattern would evaluate to True. *)
  Map[Function[r, If[Length[r] =!= 2, $Failed, Switch[Head[r],
      GreaterEqual, {"ineq", r[[1]] - r[[2]], iGreaterEqual[r[[1]], r[[2]]]},
      LessEqual, {"ineq", r[[2]] - r[[1]], iLessEqual[r[[1]], r[[2]]]},
      Equal, {"eq", r[[1]] - r[[2]], iEqual[r[[1]], r[[2]]]},
      _, $Failed]]], l]];

optimalityGen[input_, h_, spec_, vars_List, o_] := Module[{f, cons, res, m, xs, sgn, F, M,
    subs = {}, red, eqHyps = {}, hh, v, cf, sl, ineqs, gr, Fr, act, lam, G, gradF,
    sigma, st, terms, used, held, stmt, proof, rel, feas, value, lamFull},
  If[!AllTrue[vars, Head[#] === Symbol &] || vars === {}, Return[nocert[input], Module]];
  {f, cons} = If[ListQ[spec], {First[spec], optCons[And @@ Rest[spec]]}, {spec, {}}];
  If[MemberQ[cons, $Failed] || !AllTrue[Append[cons[[All, 2]], f], polyQ[#, vars] &] ||
     !SubsetQ[vars, Variables[Append[cons[[All, 2]], f]]] ||
     !FreeQ[{f, cons[[All, 2]]}, _Real | _Complex],
    Return[nocert[input], Module]];
  res = Quiet[h[spec, vars]];
  If[!MatchQ[res, {_?ratNumberQ, {__Rule}}] || !AllTrue[vars /. res[[2]], ratNumberQ],
    Return[nocert[input], Module]];
  m = res[[1]]; xs = Thread[vars -> (vars /. res[[2]])];
  sgn = If[h === Maximize, -1, 1];
  F = sgn f; M = sgn m;
  (* Eliminate linear equalities: solve for the last remaining variable with a
     constant coefficient. *)
  red = vars;
  Do[
    hh = Expand[c[[2]] /. subs];
    v = SelectFirst[Reverse[red], Function[w, Exponent[hh, w] === 1 &&
        ratNumberQ[Coefficient[hh, w, 1]] && Coefficient[hh, w, 1] =!= 0], None];
    If[v === None, If[hh === 0, Continue[], Return[nocert[input], Module]]];
    cf = Coefficient[hh, v, 1]; sl = Expand[v - hh/cf];
    subs = Append[subs /. v -> sl, v -> sl];
    AppendTo[eqHyps, {v, sl}];
    red = DeleteCases[red, v],
    {c, Select[cons, #[[1]] === "eq" &]}];
  ineqs = Select[cons, #[[1]] === "ineq" &];
  gr = Expand[ineqs[[All, 2]] /. subs];
  Fr = Expand[F /. subs];
  (* KKT: grad F = Sum lambda_j grad g_j over the active constraints at x*. *)
  act = Select[Range[Length[gr]], Expand[gr[[#]] /. xs] === 0 &];
  gradF = Expand[D[Fr, #] & /@ red] /. xs;
  lam = Which[
    red === {}, ConstantArray[0, Length[act]],
    act === {}, If[AllTrue[gradF, # === 0 &], {}, $Failed],
    True,
      G = Transpose[Table[Expand[D[gr[[j]], #] & /@ red] /. xs, {j, act}]];
      Quiet[LinearSolve[G, gradF]]];
  If[!ListQ[lam] || !AllTrue[lam, ratNumberQ] || AnyTrue[lam, # < 0 &] ||
     (act =!= {} && red =!= {} && Expand[G . lam - gradF] =!= ConstantArray[0, Length[red]]),
    Return[nocert[input], Module]];
  lamFull = ConstantArray[0, Length[gr]];
  Do[lamFull[[act[[i]]]] = lam[[i]], {i, Length[act]}];
  sigma = Expand[Fr - M - lamFull . gr];
  st = If[sigma === 0, {{}, "None"}, sosTerms[sigma, {}, red]];
  If[st === $Failed, Return[nocert[input], Module]];
  terms = st[[1]];
  (* The theorem. *)
  rel = If[h === Maximize, iLessEqual[f, m], iGreaterEqual[f, m]];
  feas = Map[# /. xs &, cons[[All, 3]] /. {(r:iGreaterEqual | iLessEqual | iEqual)[a_, b_] :>
      r[inertOf[a], inertOf[b]]}];
  value = iEqual[inertOf[f] /. xs, m];
  held = inertHeld[iAnd @@ Join[
    {If[cons === {}, iForAll[vars, rel],
      iForAll[vars, If[Length[cons] === 1, cons[[1, 3]], iAnd @@ cons[[All, 3]]], rel]]},
    feas, {value}]];
  stmt = leanThm[held, "Rat"];
  If[!StringQ[stmt], Return[nocert[input], Module]];
  used = Select[Range[Length[gr]], lamFull[[#]] =!= 0 &];
  proof = StringJoin["by",
    "\n  refine And.intro ?_ ", nestAnd[Length[feas] + 1],
    introLine[vars], If[cons === {}, "", "\n  intros"],
    Table["\n  have mth_v" <> ToString[i] <> " : " <> leanVar[eqHyps[[i, 1]]] <> " = " <>
      ratText[eqHyps[[i, 2]]] <> " := (by grind)\n  subst mth_v" <> ToString[i],
      {i, Length[eqHyps]}],
    Table["\n  have mth_g" <> ToString[j] <> " : 0 \[LessEqual] " <> ratText[gr[[j]]] <>
      " := (by grind)", {j, used}],
    sqHints[terms, {}],
    "\n  have mth_e : ", ratText[Expand[Fr - M]], " = ",
    StringRiffle[Join[
      ("(" <> ratText[lamFull[[#]]] <> ") * (" <> ratText[gr[[#]]] <> ")") & /@ used,
      If[terms === {}, {}, {sqSum[terms]}]] /. {} -> {"0"}, " + "],
    " := (by grind)", "\n  grind"];
  structGen["Optimality", input, res, held, genericStmt[stmt, "Ordered"], genericText[proof],
    <|"Value" -> m, "Minimizer" -> xs, "Multipliers" -> lamFull,
      "Constraints" -> ineqs[[All, 3]] /. $inertRules, "Eliminated" -> (Rule @@@ eqHyps),
      "Squares" -> terms, "Method" -> st[[2]], "LinearProgram" -> (sigma === 0),
      "PreCheck" -> (Expand[F /. xs] === M && AllTrue[ineqs[[All, 2]] /. xs, # >= 0 &] &&
        AllTrue[Select[cons, #[[1]] === "eq" &][[All, 2]] /. xs, # === 0 &])|>,
    {}, {}, None, $genericDomain["Ordered"]]];

(* The remaining conjuncts (feasibility, value) are closed rational facts. *)
nestAnd[1] := "(by grind)";
nestAnd[n_Integer] := "(And.intro (by grind) " <> nestAnd[n - 1] <> ")";

generate[HoldComplete[(h:Minimize | Maximize)[spec_, vs_]], o_] := With[{sv = spec},
  optimalityGen[HoldForm[h[spec, vs]], h, sv, Flatten[{vs}], o]];

(* ---------------------------------------------------------------------- *)
(* Row 15: primality -- flat Pratt certificates checked by decide +kernel   *)
(* ---------------------------------------------------------------------- *)

lucasOK[p_, a_, fs_] := PowerMod[a, p - 1, p] === 1 &&
  AllTrue[fs[[All, 1]], PowerMod[a, (p - 1)/#, p] =!= 1 &];

(* A Lucas witness for the prime p: PrimitiveRoot, else the first a >= 2. *)
prattWitness[p_, fs_] := Module[{a = Quiet[PrimitiveRoot[p]]},
  If[IntegerQ[a] && lucasOK[p, a, fs], Return[a, Module]];
  Do[If[lucasOK[p, b, fs], Return[b, Module]], {b, 2, Min[p - 1, 100000]}];
  $Failed];

(* The flat certificate of the odd prime p: steps {q, a, FactorInteger[q - 1]}
   in dependency order, deduplicated; 2 is the scheme's base. *)
prattSteps[2] = {};
prattSteps[p_Integer] := Module[{fs = FactorInteger[p - 1], a, sub},
  a = prattWitness[p, fs];
  If[a === $Failed, Return[$Failed, Module]];
  sub = prattSteps /@ DeleteCases[fs[[All, 1]], 2];
  If[MemberQ[sub, $Failed], $Failed, Append[Flatten[sub, 1], {p, a, fs}]]];
pratt[p_Integer] := With[{s = prattSteps[p]}, If[ListQ[s], DeleteDuplicates[s], $Failed]];

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

generate[HoldComplete[FactorInteger[n_]], o_] := Module[{input = HoldForm[FactorInteger[n]],
    nv = n, fi, odd, certs, pre, stmt, held},
  If[!IntegerQ[nv] || nv < 2, Return[nocert[input], Module]];
  fi = FactorInteger[nv];
  odd = DeleteCases[fi[[All, 1]], 2];
  certs = pratt /@ odd;
  If[MemberQ[certs, $Failed], Return[nocert[input], Module]];
  pre = StringJoin[Table["theorem @NAME@_pratt_" <> ToString[i] <> " : Mathilda.PrattCertifies " <>
      ToString[odd[[i]]] <> " " <> prattText[certs[[i]]] <> " := by decide +kernel\n",
    {i, Length[odd]}]];
  stmt = StringRiffle[Join[{ToString[nv] <> " = " <> StringRiffle[factorPow /@ fi, " * "]},
    ("(\[Exists] c, Mathilda.PrattCertifies " <> ToString[#] <> " c)") & /@ odd],
    " \[And] "];
  held = inertHeld[iAnd @@ Join[{iEqual[nv, iProd @@ (If[#[[2]] === 1, #[[1]], iPower @@ #] & /@ fi)]},
    iPrimeQ /@ fi[[All, 1]]]];
  Append[structGen["Primality", input, fi, held, stmt,
    "by\n  exact " <> nestAndTerms[Join[{"by decide +kernel"},
      Table["Exists.intro _ @NAME@_pratt_" <> ToString[i], {i, Length[odd]}]]],
    <|"Factorization" -> fi, "Pratt" -> AssociationThread[odd -> certs]|>,
    {$lucasStep},
    {}, {Times @@ (#[[1]]^#[[2]] & /@ fi), nv}, Integers],
    "Preamble" -> pre]];

generate[HoldComplete[PrimeQ[n_]], o_] := Module[{input = HoldForm[PrimeQ[n]], nv = n, a, cert},
  If[!IntegerQ[nv] || nv < 0, Return[nocert[input], Module]];
  Which[
    nv < 2,
      structGen["Primality", input, False, With[{x = nv}, HoldForm[!PrimeQ[x]]],
        "\[Not] Mathilda.IsPrime " <> ToString[nv], "fun h => absurd h.1 (by decide)",
        <|"Witness" -> "n < 2"|>, {}, {}, None, Integers],
    nv === 2,
      structGen["Primality", input, True, HoldForm[PrimeQ[2]], "Mathilda.IsPrime 2",
        "And.intro (by decide) (fun d h => by\n    have hle := Nat.le_of_dvd (by decide) h\n    \
rcases d with _ | _ | _ | d <;>\n      first | exact absurd (Nat.eq_zero_of_zero_dvd h) (by decide) \
| exact Or.inl rfl | exact Or.inr rfl | omega)",
        <|"Witness" -> "the divisors of 2"|>, {}, {}, None, Integers],
    PrimeQ[nv],
      cert = pratt[nv];
      If[cert === $Failed, Return[nocert[input], Module]];
      structGen["Primality", input, True, With[{x = nv}, HoldForm[PrimeQ[x]]],
        "Mathilda.PrattCertifies " <> ToString[nv] <> " " <> prattText[cert],
        "by decide +kernel", <|"Pratt" -> cert|>, {$lucasStep}, {}, None, Integers],
    True,
      a = FactorInteger[nv][[1, 1]];
      structGen["Primality", input, False, With[{x = nv}, HoldForm[!PrimeQ[x]]],
        "\[Not] Mathilda.IsPrime " <> ToString[nv],
        "fun h => absurd (h.2 " <> ToString[a] <> " (Exists.intro " <> ToString[nv/a] <>
          " (by decide +kernel))) (by decide +kernel)",
        <|"Factor" -> {a, nv/a}|>, {}, {}, None, Integers]]];

(* ---------------------------------------------------------------------- *)
(* Row 16: Sum[t, {k, a, n}] with t = polynomial(k) c^k (c rational) or a   *)
(* polynomial: Lean defines t and S by recursion from a and proves          *)
(* forall m : Nat, S m = closed(m + a) by induction.                        *)
(* ---------------------------------------------------------------------- *)

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
  If[!FreeQ[r, $Failed], Return[$Failed, Module]];
  r /. v -> If[sh === 0, ratCast[b], iPlus[ratCast[b], sh]]];

$sumBinder = Lean`Private`m1;
$termBinder = Lean`Private`k1;
$sumPass = {"Lean.Private.ratCast", "Lean.Private.m1", "Lean.Private.k1"};

sumGen[input_, t_, k_Symbol, a_, n_Symbol, cl0_, o_] := Module[{tl, cl, ttxt, ctxt, av,
    held, stmt, pre, proof, closed},
  av = a;
  (* Sum @@ {...}: a literal Sum[t, {k, ...}] here would scope (rename) the
     iterator instead of substituting the argument k. *)
  closed = If[cl0 === Automatic, If[IntegerQ[av], Sum @@ {t, {k, av, n}}, $Failed], cl0];
  (* The upper bound must be a variable, not a built-in constant: n_Symbol
     also matches Infinity, Pi, E, I, ... *)
  If[Context[n] === "System`" || MemberQ[Attributes[n], Protected | Constant],
    Return[nocert[input], Module]];
  If[!IntegerQ[av] || av < 0 || k === n || !FreeQ[closed, Sum] ||
     !SubsetQ[{k}, Union[Cases[t, _Symbol, {0, Infinity}]]] ||
     !SubsetQ[{n}, Union[Cases[closed, _Symbol, {0, Infinity}]]] ||
     !FreeQ[{t, closed}, _Real | _Complex],
    Return[nocert[input], Module]];
  (* A denominator in k (or in n in the closed form) can vanish inside the
     range, where Lean's x / 0 = 0 would give the sum a junk value; such
     sums get no certificate. *)
  If[denominatorBases[t] =!= {} && !FreeQ[denominatorBases[t], k] ||
     denominatorBases[closed] =!= {} && !FreeQ[denominatorBases[closed], n],
    Return[nocert[input], Module]];
  tl = natForm[t, k, $termBinder, 0];
  cl = natForm[closed, n, $sumBinder, av];
  If[tl === $Failed || cl === $Failed, Return[nocert[input], Module]];
  ttxt = leanStmt[inertHeld[tl], "Rat", $sumPass];
  ctxt = leanStmt[inertHeld[cl], "Rat", $sumPass];
  If[!StringQ[ttxt] || !StringQ[ctxt], Return[nocert[input], Module]];
  pre = StringJoin[
    "def @NAME@_t (k1 : Nat) : Rat := ", ttxt, "\n",
    "def @NAME@_S : Nat \[RightArrow] Rat\n",
    "  | 0 => @NAME@_t ", ToString[av], "\n",
    "  | m1 + 1 => @NAME@_S m1 + @NAME@_t (m1 + 1", If[av === 0, "", " + " <> ToString[av]],
    ")\n"];
  stmt = "\[ForAll] m1 : Nat, @NAME@_S m1 = " <> ctxt;
  proof = StringJoin["by\n  intro m1\n  induction m1 with",
    "\n  | zero => simp only [@NAME@_S, @NAME@_t] <;> (try push_cast) <;> grind",
    "\n  | succ m1 ih => simp only [@NAME@_S, @NAME@_t, ih] <;> (try push_cast) <;> grind"];
  (* Inert heads, not With: With renames the iterator of a held Sum. *)
  held = inertHeld[iForAll[n, iGreaterEqual[n, av], iEqual[iSumOp[t, {k, av, n}], closed]]];
  Append[structGen["ClosedFormSum", input, closed, held, stmt, proof,
    <|"Term" -> t, "Variable" -> k, "LowerLimit" -> av, "ClosedForm" -> closed,
      "PreCheck" -> (Expand[(closed /. n -> av) - (t /. k -> av)] === 0 &&
        Expand[(closed /. n -> n + 1) - closed - (t /. k -> n + 1)] === 0)|>,
    {"the closed form is proved for integer n >= " <> ToString[av] <> " only (Lean's \
statement is over n = " <> If[av === 0, "", ToString[av] <> " + "] <> "m with m : Nat), as \
\"Theorem\" states"},
    {}, {Sum @@ {t, {k, av, n}}, closed}],
    "Preamble" -> pre]];

(* The held input is passed whole (s, c): wrapping a held Sum in With would
   rename its iterator. *)
generate[HoldComplete[s:Sum[t_, {k_Symbol, a_, n_Symbol}]], o_] :=
  sumGen[HoldForm[s], t, k, a, n, Automatic, o];

generate[HoldComplete[c:Equal[Sum[t_, {k_Symbol, a_, n_Symbol}], rhs_]], o_] :=
  sumGen[HoldForm[c], t, k, a, n, rhs, o];


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

tier1AntiderivGen[input_, F_, f_, x_Symbol, o_] := Module[{p, ftxt, b, hb, hyps, stmt, held,
    xs, proof},
  $hdHyps = {};
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
  Join[structGen["Antiderivative", input, F, held, stmt, proof,
    <|"Form" -> "HasDerivAt", "Variable" -> x, "ProofTerm" -> p[[2]]|>, {},
    With[{d = #[[3]]}, HoldForm[d != 0]] & /@ hyps, {D[F, x], f}, Reals],
    <|"Tier" -> 1|>]];

(* Anything else: a symbol or call holding a supported shape (c = a == b;
   Certify[c]) is evaluated once and dispatched again; a fixed point that
   matches no rule has no certificate. *)
generate[HoldComplete[c_], o_] := With[{v = c},
  If[HoldComplete[v] =!= HoldComplete[c], generate[HoldComplete[v], o],
    Message[Certify::nocert, HoldForm[c]]; $Failed]];

(* ---------------------------------------------------------------------- *)
(* Certify / ProofObject                                                   *)
(* ---------------------------------------------------------------------- *)

Options[Certify] = {Engine -> "Lean", Tier -> Automatic, Domain -> Automatic,
  Method -> Automatic, TimeConstraint -> 60, Assumptions -> {}, Verify -> True,
  LeanExecutable -> Automatic, LeanProject -> Automatic};

(* HoldAll (v0.6521): the claim is held (it is certified as written), and so
   are the options until the Assumptions have been checked like the claim:
   an Evaluate or a symbol with a value in them is refused (heldOptionProblem)
   before the options are evaluated. *)
SetAttributes[Certify, HoldAll];

(* None, or $Failed after a message, for the held options of a call.  The
   Assumptions are walked held with the claim walker (sixth review, v0.6522):
   arithmetic, relations, And / Or / Not / Implies, Element[vars, domain] and
   lists are kept, and the value-preserving rational-function heads
   (Expand, Together, ...) read through, exactly as in a claim; any other
   head (Reduce, Simplify, a user function, Evaluate, D) is Certify::eval,
   and a symbol with a value (or in a reserved context) is Certify::value /
   Certify::context, before the options are evaluated -- a hypothesis
   computed by Reduce or Simplify is not the hypothesis as written.  A
   default set with SetOptions[Certify, Assumptions -> ...] was evaluated
   by SetOptions itself, so it cannot be checked as written: it is refused
   (Certify::optdefault) unless the call gives its own Assumptions. *)
SetAttributes[heldOptionProblem, HoldAllComplete];
heldOptionProblem[opts___] := Module[{as, t, bad, def},
  as = Cases[HoldComplete[opts], (Rule | RuleDelayed)[Assumptions, a_] :> HoldComplete[a],
    {1, 2}];
  If[as === {},
    def = Options[Certify, Assumptions];
    If[!MatchQ[def, {Assumptions -> ({} | True)}],
      Message[Certify::optdefault, Replace[def, {(Rule | RuleDelayed)[_, v_]} :> HoldForm[v]]];
      Return[$Failed, Module]]];
  Do[
    t = Block[{$cwAssum = True, $cwNotes = {}, $cwTrusted = {}},
      With[{w = Replace[a, HoldComplete[e_] :> cw[e]]},
        If[$cwTrusted =!= {}, cfail[D], w]]];
    bad = Cases[{t}, cvalue[s_, v_] :> {HoldForm[s], v}, Infinity];
    If[bad =!= {}, symbolMessage @@ First[bad]; Return[$Failed, Module]];
    bad = Cases[{t}, cfail[h_] :> h, Infinity];
    If[bad =!= {},
      Message[Certify::eval, Replace[a, HoldComplete[e_] :> HoldForm[Assumptions -> e]],
        First[bad]];
      Return[$Failed, Module]],
    {a, as}];
  None];

certOptions[opts_List] := Module[{keys = {Engine, Tier, Domain, Method, TimeConstraint,
    Assumptions, Verify, LeanExecutable, LeanProject}},
  AssociationThread[SymbolName /@ keys, OptionValue[Certify, opts, #] & /@ keys]];

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

(* Authentication (sixth review, v0.6522).  ProofObject[assoc] is an inert
   container: anyone can write ProofObject[<|"Verified" -> True, ...|>].  So
   every association Certify produces is registered by its Hash in this
   session, and obj["Verified"] (and the summary box) reports the stored
   verdict only for a registered object; any other ProofObject gives
   Missing["Unauthenticated"].  This guards against minted or edited objects
   in ordinary code, not against code that deliberately writes Lean`Private`
   internals; an object from another session (saved and read back) is
   unauthenticated too -- re-run Certify, or check obj["LeanSource"] with
   LeanVerify. *)
$issued[_] = False;
authenticQ[a_] := TrueQ[$issued[Hash[a]]];
verdictOf[a_] := If[authenticQ[a], Lookup[a, "Verified", Missing["UnknownProperty", "Verified"]],
  Missing["Unauthenticated"]];

(* The generator association of one item plus the tier it runs at
   ("TierUsed"), or $Failed (after a message). *)
prepare[hc_HoldComplete, o_] := Module[{g, tier, proj, rw, sides, ign},
  (* The input as written (rewriteInput); Certify::eval when it applies a
     computation Mathilda would have to evaluate first. *)
  rw = rewriteInput[hc];
  If[rw === $Failed, Return[$Failed, Module]];
  g = generate[First[rw], o];
  If[!AssociationQ[g], Return[$Failed, Module]];
  If[MemberQ[rw[[2]], Factor | FactorTerms],
    g = Join[g, <|"Kind" -> If[g["Kind"] === "Identity", "Factorization", g["Kind"]],
      "TrustedSteps" -> DeleteDuplicates[Append[g["TrustedSteps"],
        "irreducibility of the factors (not certified)"]]|>]];
  If[rw[[3]] =!= {},
    g = Append[g, "TrustedSteps" -> DeleteDuplicates[Join[g["TrustedSteps"], rw[[3]]]]]];
  (* "Input" is the input as the user wrote it (a symbol holding it shows
     its value, as before). *)
  If[!MatchQ[hc, HoldComplete[_Symbol]] && First[rw] =!= hc,
    g = Append[g, "Input" -> Replace[hc, HoldComplete[e_] :> HoldForm[e]]]];
  (* For an Equal claim, "Result" and the heuristic engine's claim are
     Mathilda's own evaluation of the two sides (what Mathilda believes); the
     Lean statement is the claim as written. *)
  If[MemberQ[{"Identity", "Factorization"}, g["Kind"]] &&
      ListQ[sides = Quiet[mathildaSides[hc]]],
    g = Join[g, <|"Claim" -> sides,
      "Result" -> With[{x = sides[[1]], y = sides[[2]]}, HoldForm[x == y]]|>]];
  (* Domain, Assumptions and Method are read by the identity generator only:
     say so when another kind ignores them (ignoring a hypothesis gives a
     stronger theorem, so it is sound, but it must not be silent). *)
  If[!MemberQ[{"Identity", "Factorization"}, g["Kind"]],
    ign = Select[{Domain -> o["Domain"], Assumptions -> o["Assumptions"],
        Method -> o["Method"]},
      !MatchQ[#[[2]], Automatic | {} | True] &];
    Do[Message[Certify::optignored, r, g["Kind"]], {r, ign}]];
  (* Tier: the generator's natural tier, unless the caller insists. *)
  tier = Switch[o["Tier"], 0, g["Tier"], 1, 1, _, g["Tier"]];
  proj = leanProject[o["LeanProject"]];
  If[tier === 1 && proj === None,
    If[o["Tier"] === 1, Message[Certify::nomathlib],
      Message[Certify::tier, g["Input"]]];
    Return[$Failed, Module]];
  If[o["Tier"] === 0 && g["Tier"] === 1,
    Message[Certify::tier, g["Input"]]; Return[$Failed, Module]];
  Append[g, "TierUsed" -> tier]];

(* Engine -> "Mathilda": PossibleZeroQ / Simplify on lhs - rhs of the claim
   (lists compare entrywise); a kind with no claim is not checked. *)
heuristicZeroQ[d_] := TrueQ[Quiet[PossibleZeroQ[d]]] || TrueQ[Quiet[Simplify[d]] === 0];
(* "Verified" is never True without Lean: the heuristic's own answer is
   carried as Missing["Heuristic", True | False]. *)
heuristicVerdict[g_] := Module[{v = Missing["NotChecked"], t = 0.},
  If[ListQ[g["Claim"]],
    t = First[AbsoluteTiming[v = Missing["Heuristic", AllTrue[
      Flatten[{g["Claim"][[1]] - g["Claim"][[2]]}], heuristicZeroQ]]]]];
  <|"Verified" -> v, "Axioms" -> {}, "Logs" -> {}, "Timing" -> t|>];

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

(* Certify a list of held items with one Lean run per tier.  The result is
   aligned with the input; an item without a certificate gives $Failed. *)
certifyItems[hcs_List, opts_List] := Module[{o, engine, gens, hb, chunks, out, exe, proj,
    idx, src, keys, todo, fresh, vs, v, lv, named},
  o = certOptions[opts];
  engine = o["Engine"];
  If[!MemberQ[{"Lean", "Mathilda"}, engine],
    Message[Certify::engine, engine]; Return[$Failed, Module]];
  (* Every Domain value has one stated meaning ($domains); anything else
     (Algebraics, "Complexes", ...) is refused rather than silently mapped. *)
  If[!AnyTrue[$domains, # === o["Domain"] &],
    Message[Certify::domain, o["Domain"]]; Return[$Failed, Module]];
  With[{why = methodProblem[o["Method"]]},
    If[why =!= None, Message[Certify::method, o["Method"], why]; Return[$Failed, Module]]];
  gens = prepare[#, o] & /@ hcs;
  hb = heartbeats[o["TimeConstraint"]];
  chunks = Table[If[AssociationQ[gens[[k]]],
      theoremChunk["mathilda_cert_" <> ToString[k], gens[[k]]["Statement"],
        gens[[k]]["Proof"], hb, Lookup[gens[[k]], "Preamble", ""]], None],
    {k, Length[gens]}];
  (* The cache keys use the generators as they are (auxiliary names still
     @NAME@); the ProofObjects show them under the theorem's actual name. *)
  named[k_] := With[{g = gens[[k]], nm = "mathilda_cert_" <> ToString[k]},
    Join[g, AssociationMap[StringReplace[g[#], "@NAME@" -> nm] &,
      Select[{"Statement", "Proof", "Preamble"}, StringQ[g[#]] &]]]];
  out = ConstantArray[$Failed, Length[gens]];
  proj = leanProject[o["LeanProject"]];
  exe = leanExe[o["LeanExecutable"]];
  Do[
    idx = Select[Range[Length[gens]],
      AssociationQ[gens[[#]]] && gens[[#]]["TierUsed"] === tier &];
    If[idx =!= {},
      src = First[assembleFile[chunks[[idx]], tier]];
      Which[
        engine === "Mathilda",
          Do[With[{g = named[k]},
              out[[k]] = mkProofObject[Append[g, "TrustedSteps" ->
                  {"Mathilda's PossibleZeroQ / Simplify (a heuristic check, not a proof)"}],
                heuristicVerdict[g], "Mathilda", "Heuristic", Missing["NotChecked"], src,
                g["Statement"]]], {k, idx}],
        !TrueQ[o["Verify"]],
          Do[out[[k]] = mkProofObject[named[k], <|"Verified" -> Missing["NotChecked"],
                "Axioms" -> {}, "Logs" -> {}, "Timing" -> 0.|>, "Lean", tier,
                Missing["NotChecked"], src, named[k]["Statement"]], {k, idx}],
        True,
          keys = AssociationThread[idx -> (cacheKey[gens[[#]], hb, exe, proj] & /@ idx)];
          todo = Select[idx, !KeyExistsQ[$certCache, keys[#]] &];
          fresh = <||>;
          If[todo =!= {},
            (* Items whose proof is a user Method script each get a Lean file
               of their own: no chunk can influence another one's verdict. *)
            With[{shared = Select[todo, !TrueQ[gens[[#]]["UserProof"]] &],
                own = Select[todo, TrueQ[gens[[#]]["UserProof"]] &]},
              If[shared =!= {},
                fresh = AssociationThread[shared -> runBatch[chunks[[shared]], tier, exe, proj]]];
              Do[fresh = Append[fresh, k -> First[runBatch[chunks[[{k}]], tier, exe, proj]]],
                {k, own}]];
            Do[If[BooleanQ[fresh[k]["Verified"]],
                $certCache = Append[$certCache, keys[k] -> fresh[k]]], {k, todo}]];
          Do[
            v = If[KeyExistsQ[fresh, k], Append[fresh[k], "Cached" -> False],
              Append[$certCache[keys[k]], "Cached" -> True]];
            lv = If[MissingQ[v["Verified"]], Missing["NotAvailable"], leanVersionOf[exe]];
            If[v["Verified"] === False, Message[Certify::rejected, gens[[k]]["Input"]]];
            out[[k]] = mkProofObject[named[k], v, "Lean", tier, lv, v["Source"],
              named[k]["Statement"]],
            {k, idx}]]],
    {tier, {0, 1}}];
  out];

(* Certify[{c1, c2, ...}] certifies a batch in one Lean run (per tier). *)
Certify[{cs___}, opts:OptionsPattern[]] := If[heldOptionProblem[opts] === $Failed, $Failed,
  certifyItems[Thread[HoldComplete[{cs}]], {opts}]];

Certify[c_, opts:OptionsPattern[]] := If[heldOptionProblem[opts] === $Failed, $Failed,
  With[{r = certifyItems[{HoldComplete[c]}, {opts}]}, If[ListQ[r], First[r], r]]];

SetAttributes[ProofObject, HoldAll];

ProofObject[e_, o:OptionsPattern[]] /; !AssociationQ[Unevaluated[e]] := Certify[e, o];

ProofObject[a_?AssociationQ]["Properties"] := Sort[Keys[a]];
ProofObject[a_?AssociationQ][k_String] :=
  If[k === "Verified", verdictOf[a], Lookup[a, k, Missing["UnknownProperty", k]]];

timingText[t_?NumericQ] := ToString[Round[t, 0.01]] <> " s";
timingText[_] := "- s";

Format[ProofObject[a_?AssociationQ]] := SequenceForm @ StringJoin[
  "ProofObject[\[LeftGuillemet]", ToString[a["Kind"]], "\[RightGuillemet] Verified: ",
  ToString[verdictOf[a]], " | ",
  If[a["Engine"] === "Lean", "Lean Tier " <> ToString[a["Tier"]], "Mathilda Heuristic"],
  " | ", timingText[a["Timing"]], "]"];

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

oracleFail[why_String] := "MATHILDA-FAIL: " <> why;

(* Strict text of an evaluated polynomial, no type ascription. *)
oracleText[e_] := With[{x = e}, toLean[HoldForm[x], None, True, {}]];

oracleVarsQ[vars_List] := AllTrue[vars, Head[#] === Symbol &];

Oracle["LinearCombination", {a_, b_}, hs_List, vars_List] /;
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

Oracle["SOS", {p_, q_}, gs_List, vars_List] /; oracleVarsQ[vars] :=
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

Oracle["Pratt", p_Integer] := Which[
  p === 2, oracleFail["2 is the base of the Pratt scheme (PrattCertifies needs an odd prime)"],
  !PrimeQ[p], oracleFail[ToString[p] <> " is not prime"],
  True, With[{c = pratt[p]}, If[ListQ[c], prattText[c],
    oracleFail["no Lucas witness found for " <> ToString[p]]]]];

Oracle["Eval", s_String] := With[{r = Quiet[ToExpression[s]]},
  If[r === $Failed, oracleFail["cannot parse " <> s], ToString[r, LeanForm]]];

Oracle[args___] := oracleFail["unsupported query " <> ToString[{args}, InputForm]];

End[];

SetAttributes[{Certify, ProofObject, LeanVerify}, Protected];

EndPackage[];
