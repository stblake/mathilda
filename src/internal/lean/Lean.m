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
PolynomialExtendedGCD, PolynomialReduce, Solve, Inverse, Det, LinearSolve, NullSpace) and \
certifies its result with the Lean 4 kernel; Certify[a == b] certifies the claim a == b and \
Certify[Implies[p1 == q1 && ..., p == q]] an ideal membership. Certify[{c1, c2, ...}] \
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
leanRename[s_String] := StringReplace[StringReplace[s, $preludeHeads],
  "\[LeftGuillemet]Lean`Private`" ~~ n:(LetterCharacter ~~ DigitCharacter..) ~~
    "\[RightGuillemet]" :> n];
leanRename[other_] := other;

(* Strict Tier 0 translation of a held statement over type ty. *)
leanStmt[held_HoldForm, ty_String] := leanRename[toLean[held, ty, True, Keys[$preludeHeads]]];

(* Lean name of a binder variable. *)
leanVar[v_Symbol] := With[{s = leanStmt[HoldForm[v], "Rat"]}, If[StringQ[s], s, "_"]];

introLine[{}] := "";
introLine[vars_List] := "\n  intro " <> StringRiffle[leanVar /@ vars, " "];

(* Inert heads keep a structured expression (p' q - p q', (q^2 + p^2) ...)
   from being normalised; inertHeld turns it into the held statement. *)
$inertRules = {iPlus -> Plus, iTimes -> Times, iPower -> Power, iEqual -> Equal,
  iUnequal -> Unequal, iAnd -> And, iOr -> Or, iImplies -> Implies,
  iEquivalent -> Equivalent};
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
  (* Quadratic surds Sqrt[d] become fresh variables s_j (relation s_j^2 = d_j). *)
  Fs = F /. $surdRules;
  bases = Union[Cases[Fs, sq[d_] :> d, {0, Infinity}]];
  svars = Table[Symbol["Lean`Private`" <>
      If[StringMatchQ[SymbolName[x], "s" ~~ DigitCharacter..], "r", "s"] <> ToString[j]],
    {j, Length[bases]}];
  sub = Thread[(sq /@ bases) -> svars];
  Fs = Fs /. sub;
  vars = Join[{x}, svars];
  If[!allowedQ[Fs, vars], Message[Certify::tier, input]; Return[$Failed, Module]];
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
  theoremChunk["mathilda_cert", g["Statement"], g["Proof"], hb]["Text"],
  g["TierUsed"], exe, If[g["TierUsed"] === 1, proj, None]}];

(* Certify a list of held items with one Lean run per tier.  The result is
   aligned with the input; an item without a certificate gives $Failed. *)
certifyItems[hcs_List, opts_List] := Module[{o, engine, gens, hb, chunks, out, exe, proj,
    idx, src, keys, todo, fresh, vs, v, lv},
  o = certOptions[opts];
  engine = o["Engine"];
  If[!MemberQ[{"Lean", "Mathilda"}, engine],
    Message[Certify::engine, engine]; Return[$Failed, Module]];
  gens = prepare[#, o] & /@ hcs;
  hb = heartbeats[o["TimeConstraint"]];
  chunks = Table[If[AssociationQ[gens[[k]]],
      theoremChunk["mathilda_cert_" <> ToString[k], gens[[k]]["Statement"],
        gens[[k]]["Proof"], hb], None], {k, Length[gens]}];
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
          Do[With[{g = gens[[k]]},
              out[[k]] = mkProofObject[Append[g, "TrustedSteps" ->
                  {"Mathilda's PossibleZeroQ / Simplify (a heuristic check, not a proof)"}],
                heuristicVerdict[g], "Mathilda", "Heuristic", Missing["NotChecked"], src,
                g["Statement"]]], {k, idx}],
        !TrueQ[o["Verify"]],
          Do[out[[k]] = mkProofObject[gens[[k]], <|"Verified" -> Missing["NotChecked"],
                "Axioms" -> {}, "Logs" -> {}, "Timing" -> 0.|>, "Lean", tier,
                Missing["NotChecked"], src, gens[[k]]["Statement"]], {k, idx}],
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
            out[[k]] = mkProofObject[gens[[k]], v, "Lean", tier, lv, v["Source"],
              gens[[k]]["Statement"]],
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
