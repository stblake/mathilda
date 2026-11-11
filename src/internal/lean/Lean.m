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
(* A generator may also supply a "Preamble" (auxiliary defs / lemmas placed
   before the theorem).  The statement, proof and preamble name their
   auxiliary declarations @NAME@_..., replaced here by the theorem's name. *)
theoremChunk[name_String, stmt_String, proof_String, hb_Integer, pre_String:""] := <|
  "Name" -> name,
  "Text" -> StringReplace[StringJoin[pre,
    If[hb > 0, "set_option maxHeartbeats " <> ToString[hb] <> " in\n", ""],
    "theorem ", name, " :\n    ", stmt, " := ", proof, "\n",
    "#print axioms ", name, "\n"], "@NAME@" -> name]|>;

fileHeader[tier_] := StringJoin[
  If[tier === 1, "import Mathlib\nset_option linter.unusedTactic false\n\
set_option linter.unreachableTactic false\n\n", ""],
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

(* ---------------------------------------------------------------------- *)
(* Shared helpers for the structured certificates (rows 4-10)              *)
(* ---------------------------------------------------------------------- *)

(* Prelude functions in statements are written with private placeholder
   heads and renamed after translation.  A symbol in a context on
   $ContextPath (Mathilda`matMul) would capture the user's short name. *)
$preludeHeads = {"Lean.Private.matMul" -> "Mathilda.matMul",
  "Lean.Private.matVec" -> "Mathilda.matVec", "Lean.Private.idMat" -> "Mathilda.idMat",
  "Lean.Private.det" -> "Mathilda.det"};

(* Private generator variables (the surd symbols s1, s2, ...) print as
   <<Lean`Private`s1>>; they are renamed to their bare names. *)
leanRename[s_String] := StringReplace[StringReplace[StringReplace[s, $preludeHeads],
  "\[LeftGuillemet]Lean`Private`" ~~ n:(LetterCharacter ~~ DigitCharacter..) ~~
    "\[RightGuillemet]" :> n],
  (* Lean`Private`ratCast[m1]: the Nat binder m1 used as a rational value. *)
  "Lean.Private.ratCast " ~~ n:(LetterCharacter ~~ DigitCharacter..) :>
    "(" <> n <> " : Rat)"];
leanRename[other_] := other;

(* Strict Tier 0 translation of a held statement over type ty. *)
leanStmt[held_HoldForm, ty_String, extra_List:{}] :=
  leanRename[toLean[held, ty, True, Join[Keys[$preludeHeads], extra]]];

(* Lean name of a binder variable. *)
leanVar[v_Symbol] := With[{s = leanStmt[HoldForm[v], "Rat"]}, If[StringQ[s], s, "_"]];

introLine[{}] := "";
introLine[vars_List] := "\n  intro " <> StringRiffle[leanVar /@ vars, " "];

(* Inert heads keep a structured expression (p' q - p q', (q^2 + p^2) ...)
   from being normalised; inertHeld turns it into the held statement. *)
$inertRules = {iPlus -> Plus, iTimes -> Times, iPower -> Power, iEqual -> Equal,
  iUnequal -> Unequal, iAnd -> And, iOr -> Or, iImplies -> Implies,
  iEquivalent -> Equivalent, iGreaterEqual -> GreaterEqual, iLessEqual -> LessEqual,
  iNot -> Not, iForAll -> ForAll, iPrimeQ -> PrimeQ, iSumOp -> Sum};
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

antiderivGen[input_, f_, x_Symbol, o_] := Module[{F, tf, a, b, Fs, bases, svars, vars,
    sub, terms, rpart = 0, elems = {}, pieces, tr, r, w, dens, lhs, lhsV, rels, red,
    ks, rem, rhs, held, stmt, nonzero, hyps, back, fs, ph, cf},
  F = If[$antiderivOverride =!= None, $antiderivOverride, Integrate[f, x]];
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
  svars = Table[Symbol["Lean`Private`" <>
      If[StringMatchQ[SymbolName[x], "s" ~~ DigitCharacter..], "r", "s"] <> ToString[j]],
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
  stmt = leanStmt[held, "Rat"];
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
     "Together normal form of the integrand (f = a/b) and of each argument u_i = p_i/q_i",
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
  stmt = leanStmt[held, ty];
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
  stmt = leanStmt[held, ty];
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
  ty = autoDomain[{ps, q, cs}];
  stmt = leanStmt[held, ty];
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
  structGen["IdealMembership", input, cs, held, stmt, proof,
    <|"Cofactors" -> cs, "Generators" -> ps, "Goal" -> q|>, {}, {}, {q, cs . ps},
    leanType[ty]]];

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

solveGen[input_, eqs_List, vars_List, sol_, o_] := Module[{x, p, roots, lc, mult, fact,
    held, stmt, ftxt, etxt, l, r, alts, sound},
  If[!MatchQ[sol, {{__Rule} ..}] || !AllTrue[Flatten[vars /. sol], ratNumberQ],
    Return[nocert[input], Module]];
  (* Completeness: one univariate polynomial equation that splits over Q. *)
  If[Length[vars] === 1 && Length[eqs] === 1,
    x = First[vars]; {l, r} = First[eqs]; p = Expand[l - r];
    If[polyQ[p, {x}] && SubsetQ[{x}, Variables[p]] && Exponent[p, x] > 0,
      roots = DeleteDuplicates[x /. sol];
      lc = Coefficient[p, x, Exponent[p, x]];
      mult = rootMultiplicity[p, x, #] & /@ roots;
      fact = iProd[lc, Sequence @@ MapThread[
        If[#2 === 1, iPlus[x, -#1], iPower[iPlus[x, -#1], #2]] &, {roots, mult}]];
      If[Expand[p - inertValue[fact]] === 0,
        alts = iEqual[x, #] & /@ roots;
        held = forAllHeld[inertHeld[iEquivalent[iEqual[l, r],
            If[Length[alts] === 1, First[alts], iOr @@ alts]]], {x}];
        stmt = leanStmt[held, "Rat"];
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
  (* Soundness: substitute every solution into every equation (held, so the
     substitution reaches Lean) and let the kernel evaluate. *)
  sound = With[{e = If[Length[eqs] === 1, iEqual @@ First[eqs], iAnd @@ (iEqual @@@ eqs)]},
    HoldForm[e]];
  held = inertHeld[If[Length[sol] === 1, sound /. First[sol], iAnd @@ ((sound /. #) & /@ sol)]];
  stmt = leanStmt[held, "Rat"];
  If[!StringQ[stmt], Return[nocert[input], Module]];
  structGen["Roots", input, sol, held, stmt, "by\n  decide +kernel",
    <|"Completeness" -> False, "Solutions" -> sol|>,
    {"completeness: that Solve returned every solution (only soundness -- each solution \
satisfies each equation -- is certified)"},
    {}, {Flatten[Table[(#[[1]] - #[[2]]) & /@ eqs /. s, {s, sol}]], 0}]];

generate[HoldComplete[Solve[eqs_, vs_]], o_] := Module[{input = HoldForm[Solve[eqs, vs]],
    el, vars, sol},
  el = Cases[HoldComplete[eqs] /. {HoldPattern[And[e___]] :> {e}},
    HoldPattern[Equal[a_, b_]] :> {a, b}, {0, 2}];
  vars = Flatten[{vs}];
  If[el === {} || !AllTrue[vars, Head[#] === Symbol &], Return[nocert[input], Module]];
  sol = Quiet[Solve[eqs, vs]];
  solveGen[input, el, vars, sol, o]];

(* ---------------------------------------------------------------------- *)
(* Rows 8-10: Inverse, Det (n <= 6), LinearSolve, NullSpace over Q,        *)
(* checked by the kernel evaluating the prelude's matMul / det / matVec    *)
(* ---------------------------------------------------------------------- *)

matrixGen[kind_, input_, result_, held_, cert_, trusted_, claim_] := With[
  {stmt = leanStmt[held, "Rat"]},
  If[!StringQ[stmt], nocert[input],
    structGen[kind, input, result, held, stmt, "by\n  decide +kernel", cert, trusted, {},
      claim]]];

generate[HoldComplete[Inverse[m_]], o_] := Module[{input = HoldForm[Inverse[m]], a = m,
    b, n},
  If[!ratMatrixQ[a] || Length[a] =!= Length[First[a]], Return[nocert[input], Module]];
  b = Quiet[Inverse[a]];
  If[!ratMatrixQ[b], Return[nocert[input], Module]];
  n = Length[a];
  matrixGen["MatrixInverse", input, b,
    With[{aa = a, bb = b, k = n},
      HoldForm[matMul[aa, bb] == idMat[k] && matMul[bb, aa] == idMat[k]]],
    <|"Matrix" -> a, "Inverse" -> b|>, {}, {{a . b, b . a}, {IdentityMatrix[n], IdentityMatrix[n]}}]];

generate[HoldComplete[Det[m_]], o_] := Module[{input = HoldForm[Det[m]], a = m, d},
  If[!ratMatrixQ[a] || Length[a] =!= Length[First[a]] || Length[a] > 6,
    Return[nocert[input], Module]];
  d = Det[a];
  matrixGen["Determinant", input, d, With[{aa = a, dd = d}, HoldForm[det[aa] == dd]],
    <|"Matrix" -> a, "Method" -> "Laplace expansion along row 0 (prelude Mathilda.det)"|>,
    {}, {d, Det[a]}]];

generate[HoldComplete[LinearSolve[m_, v_]], o_] := Module[{input = HoldForm[LinearSolve[m, v]],
    a = m, b = v, x},
  If[!ratMatrixQ[a] || !ratVectorQ[b] || Length[b] =!= Length[a],
    Return[nocert[input], Module]];
  x = Quiet[LinearSolve[a, b]];
  If[!ratVectorQ[x], Return[nocert[input], Module]];
  matrixGen["LinearSolve", input, x, With[{aa = a, xx = x, bb = b}, HoldForm[matVec[aa, xx] == bb]],
    <|"Matrix" -> a, "RightHandSide" -> b, "Solution" -> x|>, {}, {a . x, b}]];

generate[HoldComplete[NullSpace[m_]], o_] := Module[{input = HoldForm[NullSpace[m]], a = m,
    vs, zero},
  If[!ratMatrixQ[a], Return[nocert[input], Module]];
  vs = NullSpace[a];
  If[!MatchQ[vs, {__}] || !AllTrue[vs, ratVectorQ], Return[nocert[input], Module]];
  zero = ConstantArray[0, Length[a]];
  matrixGen["NullSpace", input, vs,
    inertHeld[If[Length[vs] === 1, First[#], iAnd @@ #] &[
      With[{aa = a, z = zero}, iEqual[matVec[aa, #], z]] & /@ vs]],
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

$overQSteps = {"the statement is checked over Q (Tier 0); the certificate identity proved \
inside it is a polynomial identity, so the bound also holds over R (that transfer is not \
checked by Lean)"};

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
  stmt = leanStmt[held, "Rat"];
  If[!StringQ[stmt], Return[nocert[input], Module]];
  ghints = Table["\n  have mth_g" <> ToString[i] <> " : 0 \[LessEqual] " <> ratText[used[[i]]] <>
    " := (by grind)", {i, Length[used]}];
  proof = StringJoin["by", introLine[vars], If[conds === {}, "", "\n  intros"], ghints,
    sqHints[terms, gnames],
    "\n  have mth_e : ", ratText[target], " = ", If[terms === {}, "0", sqSum[terms]],
    " := (by grind)", "\n  grind"];
  structGen["Nonnegativity", input, True, held, stmt, proof,
    <|"Squares" -> terms, "Method" -> st[[2]], "Target" -> target, "Multipliers" -> used,
      "PreCheck" -> True|>,
    $overQSteps, {}, None]];

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
  stmt = leanStmt[held, "Rat"];
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
  structGen["Optimality", input, res, held, stmt, proof,
    <|"Value" -> m, "Minimizer" -> xs, "Multipliers" -> lamFull,
      "Constraints" -> ineqs[[All, 3]] /. $inertRules, "Eliminated" -> (Rule @@@ eqHyps),
      "Squares" -> terms, "Method" -> st[[2]], "LinearProgram" -> (sigma === 0),
      "PreCheck" -> (Expand[F /. xs] === M && AllTrue[ineqs[[All, 2]] /. xs, # >= 0 &] &&
        AllTrue[Select[cons, #[[1]] === "eq" &][[All, 2]] /. xs, # === 0 &])|>,
    $overQSteps, {}, None]];

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
        ToString[nv] <> " < 2", "by decide", <|"Witness" -> "n < 2"|>, {}, {}, None, Integers],
    nv === 2,
      structGen["Primality", input, True, HoldForm[PrimeQ[2]],
        "\[ForAll] d, d < 3 \[RightArrow] Dvd.dvd d 2 \[RightArrow] d = 1 \[Or] d = 2",
        "by decide", <|"Witness" -> "divisors of 2 below 3"|>,
        {"a divisor of n is at most n (so checking d < 3 suffices)"}, {}, None, Integers],
    PrimeQ[nv],
      cert = pratt[nv];
      If[cert === $Failed, Return[nocert[input], Module]];
      structGen["Primality", input, True, With[{x = nv}, HoldForm[PrimeQ[x]]],
        "Mathilda.PrattCertifies " <> ToString[nv] <> " " <> prattText[cert],
        "by decide +kernel", <|"Pratt" -> cert|>, {$lucasStep}, {}, None, Integers],
    True,
      a = FactorInteger[nv][[1, 1]];
      structGen["Primality", input, False, With[{x = nv}, HoldForm[!PrimeQ[x]]],
        ToString[nv] <> " = " <> ToString[a] <> " * " <> ToString[nv/a] <> " \[And] 1 < " <>
          ToString[a] <> " \[And] " <> ToString[a] <> " < " <> ToString[nv],
        "by decide +kernel", <|"Factor" -> {a, nv/a}|>, {}, {}, None, Integers]]];

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
  If[!IntegerQ[av] || av < 0 || k === n || !FreeQ[closed, Sum] ||
     !SubsetQ[{k}, Union[Cases[t, _Symbol, {0, Infinity}]]] ||
     !SubsetQ[{n}, Union[Cases[closed, _Symbol, {0, Infinity}]]] ||
     !FreeQ[{t, closed}, _Real | _Complex],
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
    {}, {}, {Sum @@ {t, {k, av, n}}, closed}],
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
realConst[c_] := With[{s = toLean[HoldForm[c], "Real", False, {}]},
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
  "Timing" -> v["Timing"], "Cached" -> TrueQ[v["Cached"]], "LeanVersion" -> lv|>},
  ProofObject[a]];  (* With: the held slot receives the evaluated Association *)

(* The generator association of one item plus the tier it runs at
   ("TierUsed"), or $Failed (after a message). *)
prepare[hc_HoldComplete, o_] := Module[{g, tier, proj},
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
  Append[g, "TierUsed" -> tier]];

(* Engine -> "Mathilda": PossibleZeroQ / Simplify on lhs - rhs of the claim
   (lists compare entrywise); a kind with no claim is not checked. *)
heuristicZeroQ[d_] := TrueQ[Quiet[PossibleZeroQ[d]]] || TrueQ[Quiet[Simplify[d]] === 0];
heuristicVerdict[g_] := Module[{v = Missing["NotChecked"], t = 0.},
  If[ListQ[g["Claim"]],
    t = First[AbsoluteTiming[v = AllTrue[
      Flatten[{g["Claim"][[1]] - g["Claim"][[2]]}], heuristicZeroQ]]]];
  <|"Verified" -> v, "Axioms" -> {}, "Logs" -> {}, "Timing" -> t|>];

(* Session cache of Lean verdicts, keyed by the theorem text under a fixed
   name (so a theorem is found whatever its position in a batch), the tier,
   the executable and the project.  Only definitive verdicts are stored. *)
$certCache = <||>;
cacheKey[g_, hb_, exe_, proj_] := Hash[{
  theoremChunk["mathilda_cert", g["Statement"], g["Proof"], hb,
    Lookup[g, "Preamble", ""]]["Text"],
  g["TierUsed"], exe, If[g["TierUsed"] === 1, proj, None]}];

(* Certify a list of held items with one Lean run per tier.  The result is
   aligned with the input; an item without a certificate gives $Failed. *)
certifyItems[hcs_List, opts_List] := Module[{o, engine, gens, hb, chunks, out, exe, proj,
    idx, src, keys, todo, fresh, vs, v, lv, named},
  o = certOptions[opts];
  engine = o["Engine"];
  If[!MemberQ[{"Lean", "Mathilda"}, engine],
    Message[Certify::engine, engine]; Return[$Failed, Module]];
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
            vs = runBatch[chunks[[todo]], tier, exe, proj];
            fresh = AssociationThread[todo -> vs];
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
Certify[{cs___}, opts:OptionsPattern[]] := certifyItems[Thread[HoldComplete[{cs}]], {opts}];

Certify[c_, opts:OptionsPattern[]] := With[{r = certifyItems[{HoldComplete[c]}, {opts}]},
  If[ListQ[r], First[r], r]];

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
