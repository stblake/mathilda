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
Together, Apart, Cancel, Simplify, FullSimplify) and certifies its result with the Lean 4 \
kernel; Certify[a == b] certifies the claim a == b. The result is a ProofObject whose \
\"Verified\" property is True only when Lean accepted the theorem.";

ProofObject::usage =
  "ProofObject[<|...|>] is a Lean 4 certificate produced by Certify; obj[\"prop\"] gives a \
property and obj[\"Properties\"] lists them. ProofObject[comp, opts] (comp not an \
Association) is the same as Certify[comp, opts].";

LeanVerify::usage =
  "LeanVerify[src] checks the Lean 4 source string src with the Lean executable and gives \
<|\"Success\", \"ExitCode\", \"Messages\", \"Axioms\", \"Timing\", \"Command\"|>. Success \
requires exit code 0, no error, no sorry, and only the standard axioms for every theorem.";

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
\"Mathilda\" (a heuristic PossibleZeroQ / Simplify check, for contrast).";

Tier::usage =
  "Tier is an option of Certify: 0 (core Lean), 1 (Mathlib, needs $LeanProject) or \
Automatic.";

Domain::usage =
  "Domain is an option of Certify: Integers, Rationals, Reals (Tier 1) or Automatic.";

Verify::usage =
  "Verify is an option of Certify: False generates the Lean source without running Lean \
(\"Verified\" -> Missing[\"NotChecked\"]).";

LeanExecutable::usage =
  "LeanExecutable is an option of Certify and LeanVerify that overrides $LeanExecutable.";

LeanProject::usage =
  "LeanProject is an option of Certify and LeanVerify that overrides $LeanProject.";

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

(* Names declared with `theorem` in src. *)
declaredTheorems[src_String] := DeleteDuplicates[StringCases[src,
   RegularExpression["(?m)^\\s*(?:private\\s+|protected\\s+)?theorem\\s+([^\\s:(){}\\[\\]]+)"] ->
    "$1"]];

(* src with a `#print axioms` line for every declared theorem that lacks one. *)
withAxiomReports[src_String] := Module[{missing},
  missing = Select[declaredTheorems[src],
    !StringContainsQ[src, "#print axioms " <> #] &];
  If[missing === {}, src,
    StringJoin[src, If[StringEndsQ[src, "\n"], "", "\n"],
      StringJoin[("#print axioms " <> # <> "\n") & /@ missing]]]];

Options[LeanVerify] = {LeanExecutable -> Automatic, LeanProject -> Automatic,
  KeepFiles -> False};

LeanVerify[src_String, opts:OptionsPattern[]] := Module[{exe, proj, full, r, names, ok},
  exe = leanExe[OptionValue[LeanExecutable]];
  proj = leanProject[OptionValue[LeanProject]];
  full = withAxiomReports[src];
  r = Lean`Internal`VerifySource[full, exe, proj, TrueQ[OptionValue[KeepFiles]]];
  If[!AssociationQ[r], Return[$Failed, Module]];
  If[missingLeanQ[r], Message[LeanVerify::nolean]];
  names = declaredTheorems[src];
  ok = TrueQ[r["Success"]] &&
    AllTrue[names, ListQ[Lookup[r["Axioms"], #, None]] &];
  Join[<|"Success" -> ok|>, KeyDrop[r, "Success"]]];

(* ---------------------------------------------------------------------- *)
(* Translation and source assembly                                         *)
(* ---------------------------------------------------------------------- *)

leanType["Int"] = Integers; leanType["Rat"] = Rationals; leanType["Real"] = Reals;

(* Strict translation of a held statement; $Failed (quietly) when refused. *)
toLean[held_HoldForm, ty_, coreOnly_, pass_List:{}] :=
  Quiet[Lean`Internal`ToLean[held, ty, coreOnly, pass]];

heartbeats[t_?NumericQ] /; t > 0 := Round[40000 t];
heartbeats[_] := 0;

proofScript[method_] := Which[
  method === Automatic || method === "grind", "by\n  intros\n  grind",
  method === "decide", "by\n  decide",
  method === "omega", "by\n  intros\n  omega",
  StringQ[method], "by\n  " <> method,
  True, "by\n  intros\n  grind"];

(* One theorem chunk: <|"Name", "Text"|>.  The #print axioms report sits in
   the chunk so an unparsable theorem's report error is attributed to it. *)
theoremChunk[name_String, stmt_String, proof_String, hb_Integer] := <|
  "Name" -> name,
  "Text" -> StringJoin[
    If[hb > 0, "set_option maxHeartbeats " <> ToString[hb] <> " in\n", ""],
    "theorem ", name, " :\n    ", stmt, " := ", proof, "\n",
    "#print axioms ", name, "\n"]|>;

fileHeader[tier_] := StringJoin[
  If[tier === 1, "import Mathlib\n\n", ""],
  "/- Mathilda ", ToString[$VersionNumber], " certificate (",
  If[tier === 1, "Tier 1, Mathlib", "Tier 0, core Lean"],
  "). Generated: Lean is the judge; Mathilda is an untrusted oracle. -/\n\n"];

(* The file for a list of chunks plus each chunk's line span {first, last}. *)
assembleFile[chunks_List, tier_] := Module[{head, line, spans = {}, body = {}},
  head = fileHeader[tier] <> prelude[] <> "\n";
  line = StringCount[head, "\n"] + 1;
  Do[
    With[{t = c["Text"]},
      AppendTo[spans, {line, line + StringCount[t, "\n"] - 1}];
      AppendTo[body, t <> "\n"];
      line += StringCount[t, "\n"] + 1],
    {c, chunks}];
  {StringJoin[head, body], spans}];

(* ---------------------------------------------------------------------- *)
(* Running Lean and attributing verdicts (the Verified contract)           *)
(* ---------------------------------------------------------------------- *)

badMessageQ[m_] := m["Severity"] === "error" || m["Kind"] === "hasSorry";
inSpanQ[m_, {a_, b_}] := a <= m["Line"] <= b;

(* Verdict of the theorem `name` with line span `span`, given the spans of all
   theorems in the file.  True needs: no error / sorry in its span, none
   outside every span (prelude, header), an axiom report naming only the
   standard axioms, and -- when the exit code is non-zero -- every error
   attributable to some theorem. *)
attribute[r_, name_, span_, spans_] := Module[{bad, own, stray, ax, logs},
  bad = Select[r["Messages"], badMessageQ];
  own = Select[bad, inSpanQ[#, span] &];
  stray = Select[bad, Function[m, NoneTrue[spans, inSpanQ[m, #] &]]];
  ax = Lookup[r["Axioms"], name, Missing["NotReported"]];
  logs = Function[m, <|"Severity" -> m["Severity"], "Line" -> m["Line"],
      "Text" -> m["Text"]|>] /@
    Select[r["Messages"], Function[m, inSpanQ[m, span] ||
        (badMessageQ[m] && NoneTrue[spans, inSpanQ[m, #] &])]];
  <|"Verified" -> (own === {} && stray === {} && ListQ[ax] &&
      SubsetQ[$standardAxioms, ax] && (r["ExitCode"] === 0 || bad =!= {})),
    "Axioms" -> If[ListQ[ax], ax, {}], "Logs" -> logs|>];

(* Run one file of chunks; a list of per-chunk verdict associations. *)
runBatch[chunks_List, tier_, exe_, proj_] := Module[{src, spans, r},
  {src, spans} = assembleFile[chunks, tier];
  r = Lean`Internal`VerifySource[src, exe, If[tier === 1, proj, None], False];
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
      MapThread[Join[attribute[r, #1["Name"], #2, spans],
          <|"Timing" -> r["Timing"], "Source" -> src|>] &, {chunks, spans}]]];

(* ---------------------------------------------------------------------- *)
(* Certificate generators (no Lean call)                                   *)
(* ---------------------------------------------------------------------- *)

(* Domain selection: Integers when the statement has integer coefficients and
   no division, Rationals otherwise. *)
autoDomain[stmt_] := If[FreeQ[stmt, _Rational] && FreeQ[stmt, Power[_, _?Negative]],
  "Int", "Rat"];

domainType[Automatic, stmt_] := autoDomain[stmt];
domainType[Integers, stmt_] := "Int";
domainType[Rationals, stmt_] := "Rat";
domainType[Reals, stmt_] := "Real";
domainType[_, stmt_] := autoDomain[stmt];

polyQ[e_, vars_] := PolynomialQ[e, vars] && FreeQ[e, Power[_, _?Negative]];

(* The certified statement for the identity a == b over `vars`:
   {held statement, hypotheses, form, extra certificate data, extra trusted steps}. *)
identityStatement[a_, b_, vars_] := Module[{ta, tb, na, da, nb, db},
  If[polyQ[a, vars] && polyQ[b, vars],
    Return[{HoldForm[a == b], {}, "Polynomial", <||>, {}}, Module]];
  ta = Together[a]; tb = Together[b];
  {na, da} = {Numerator[ta], Denominator[ta]};
  {nb, db} = {Numerator[tb], Denominator[tb]};
  If[!AllTrue[{na, da, nb, db}, polyQ[#, vars] &], Return[$Failed, Module]];
  With[{l = If[db === 1, na, na db], r = If[da === 1, nb, nb da]},
    {HoldForm[l == r],
     Join[If[da === 1, {}, With[{d = da}, {HoldForm[d != 0]}]],
          If[db === 1, {}, With[{d = db}, {HoldForm[d != 0]}]]],
     "ClearedDenominators",
     <|"Numerators" -> {na, nb}, "Denominators" -> {da, db}|>,
     {"Together normal form of each side (num/den used to clear denominators)"}}]];

(* Generator for an identity a == b (rows 1-2).  Returns the generator
   association, or $Failed after Certify::nocert. *)
identityGen[kind_, input_HoldForm, result_, a_, b_, trusted_List, o_] :=
  Module[{vars, st, ty, held, assum, lean, tier},
  vars = Sort[Variables[{a, b}]];
  If[!AllTrue[vars, Head[#] === Symbol &],
    Message[Certify::nocert, input]; Return[$Failed, Module]];
  st = identityStatement[a, b, vars];
  If[st === $Failed, Message[Certify::nocert, input]; Return[$Failed, Module]];
  ty = domainType[o["Domain"], First[st]];
  tier = If[ty === "Real", 1, 0];
  assum = Flatten[{o["Assumptions"]}];
  (* Bind the variables (and assumptions) around the held statement without
     ever evaluating it: Equal of two identical sides would become True. *)
  held = With[{vs = vars, h = If[Length[assum] === 1, First[assum], And @@ assum]},
    Which[
      vs === {} && assum === {}, First[st],
      vs === {}, First[st] /. HoldForm[s_] :> HoldForm[Implies[h, s]],
      assum === {}, First[st] /. HoldForm[s_] :> HoldForm[ForAll[vs, s]],
      True, First[st] /. HoldForm[s_] :> HoldForm[ForAll[vs, h, s]]]];
  lean = toLean[held, ty, tier === 0];
  If[!StringQ[lean], Message[Certify::nocert, input]; Return[$Failed, Module]];
  <|"Kind" -> kind, "Input" -> input, "Result" -> result, "Theorem" -> held,
    "Statement" -> lean, "Proof" -> proofScript[o["Method"]],
    "Certificate" -> Join[<|"Form" -> st[[3]], "Variables" -> vars|>, st[[4]]],
    "TrustedSteps" -> Join[trusted, st[[5]]], "Hypotheses" -> st[[2]],
    "Tier" -> tier, "Domain" -> leanType[ty], "Claim" -> {a, b}|>];

$identityHeads = {Expand, Collect, Together, Apart, Cancel, Simplify, FullSimplify};
$factorHeads = {Factor, FactorTerms};

(* Dispatcher: one rule per supported computation / claim shape. *)
generate[HoldComplete[Equal[a_, b_]], o_] := Module[{av = a, bv = b},
  With[{x = av, y = bv},
    identityGen["Identity", HoldForm[a == b], HoldForm[x == y], av, bv, {}, o]]];

generate[HoldComplete[(h_Symbol)[e_, rest___]], o_] /;
    MemberQ[Join[$identityHeads, $factorHeads], h] :=
  Module[{ev = e, res},
    res = h[ev, rest];
    identityGen[
      If[MemberQ[$factorHeads, h], "Factorization", "Identity"],
      HoldForm[h[e, rest]], res, ev, res,
      If[MemberQ[$factorHeads, h], {"irreducibility of the factors (not certified)"}, {}],
      o]];

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

SetAttributes[Certify, HoldFirst];

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
  "Timing" -> v["Timing"], "Cached" -> False, "LeanVersion" -> lv|>},
  ProofObject[a]];  (* With: the held slot receives the evaluated Association *)

certify[hc_HoldComplete, opts_List] := Module[{o, g, engine, tier, chunk, exe, proj,
    v, src, lv, t},
  o = certOptions[opts];
  engine = o["Engine"];
  If[!MemberQ[{"Lean", "Mathilda"}, engine],
    Message[Certify::engine, engine]; Return[$Failed, Module]];
  g = generate[hc, o];
  If[!AssociationQ[g], Return[$Failed, Module]];
  (* Tier: the generator's natural tier, unless the caller insists. *)
  tier = Switch[o["Tier"], 0, g["Tier"], 1, 1, _, g["Tier"]];
  proj = leanProject[o["LeanProject"]];
  If[tier === 1 && proj === None,
    If[o["Tier"] === 1, Message[Certify::nomathlib],
      Message[Certify::tier, g["Input"]]];
    Return[$Failed, Module]];
  If[o["Tier"] === 0 && g["Tier"] === 1,
    Message[Certify::tier, g["Input"]]; Return[$Failed, Module]];
  chunk = theoremChunk["mathilda_cert_1", g["Statement"], g["Proof"],
    heartbeats[o["TimeConstraint"]]];
  src = First[assembleFile[{chunk}, tier]];
  Which[
    engine === "Mathilda",
      t = First[AbsoluteTiming[
        v = With[{d = g["Claim"][[1]] - g["Claim"][[2]]},
          TrueQ[Quiet[PossibleZeroQ[d]]] || TrueQ[Quiet[Simplify[d]] === 0]]]];
      With[{gh = Append[g, "TrustedSteps" ->
            {"Mathilda's PossibleZeroQ / Simplify (a heuristic check, not a proof)"}]},
        mkProofObject[gh, <|"Verified" -> v, "Axioms" -> {}, "Logs" -> {}, "Timing" -> t|>,
          "Mathilda", "Heuristic", Missing["NotChecked"], src, g["Statement"]]],
    !TrueQ[o["Verify"]],
      mkProofObject[g, <|"Verified" -> Missing["NotChecked"], "Axioms" -> {}, "Logs" -> {},
          "Timing" -> 0.|>, "Lean", tier, Missing["NotChecked"], src, g["Statement"]],
    True,
      exe = leanExe[o["LeanExecutable"]];
      v = First[runBatch[{chunk}, tier, exe, proj]];
      lv = If[MissingQ[v["Verified"]], Missing["NotAvailable"], leanVersionOf[exe]];
      If[v["Verified"] === False, Message[Certify::rejected, g["Input"]]];
      mkProofObject[g, v, "Lean", tier, lv, v["Source"], g["Statement"]]]];

Certify[c_, opts:OptionsPattern[]] := certify[HoldComplete[c], {opts}];

SetAttributes[ProofObject, HoldFirst];

ProofObject[e_, o:OptionsPattern[]] /; !AssociationQ[Unevaluated[e]] := Certify[e, o];

ProofObject[a_?AssociationQ]["Properties"] := Sort[Keys[a]];
ProofObject[a_?AssociationQ][k_String] := Lookup[a, k, Missing["UnknownProperty", k]];

timingText[t_?NumericQ] := ToString[Round[t, 0.01]] <> " s";
timingText[_] := "- s";

Format[ProofObject[a_?AssociationQ]] := SequenceForm @ StringJoin[
  "ProofObject[\[LeftGuillemet]", ToString[a["Kind"]], "\[RightGuillemet] Verified: ",
  ToString[a["Verified"]], " | ",
  If[a["Engine"] === "Lean", "Lean Tier " <> ToString[a["Tier"]], "Mathilda Heuristic"],
  " | ", timingText[a["Timing"]], "]"];

End[];

SetAttributes[{Certify, ProofObject, LeanVerify}, Protected];

EndPackage[];
