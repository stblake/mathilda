(* ====================================================================== *)
(* Lean.m -- Certify, ProofObject, LeanVerify: the prover-neutral core.   *)
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
(* This file: claim walker, domains, generators (neutral IR), backend      *)
(* registry, pipeline, ProofObject, heuristic engine.  LeanBackend.m (loaded *)
(* at the end): the Lean backend (discovery, LeanVerify, encoders, runs).   *)
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
is True only when Lean accepted the theorem. Options (see Options[Certify]; a name is matched \
by its text, so \"Engine\" -> v and Engine -> v are the same): \"Engine\" (\"Lean\" or \
\"Mathilda\", a heuristic check whose \"Verified\" is never True), \"Tier\" (0, 1 or \
Automatic), \"Domain\" (Automatic, Complexes, Reals, Rationals, Integers), Method (a tactic \
script), TimeConstraint, Assumptions, \"Verify\" (False: generate the source only), \
\"LeanExecutable\", \"LeanProject\".";

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
sorry, ...) is refused unrun. Options: \"LeanExecutable\", \"LeanProject\" (override \
$LeanExecutable / $LeanProject) and \"KeepFiles\" (True keeps the .lean file, returned as \
\"SourceFile\").";

$LeanExecutable::usage =
  "$LeanExecutable is the Lean 4 executable used by Certify and LeanVerify: Automatic \
(MATHILDA_LEAN, else the newest ~/.elan toolchain, else lean on PATH) or a path.";

$LeanProject::usage =
  "$LeanProject is None or the directory of a Lake project with Mathlib built, used for \
Tier -> 1 certificates. It defaults to the environment variable MATHILDA_LEAN_PROJECT.";

$LeanVersion::usage =
  "$LeanVersion gives the version string of the Lean executable, or \
Missing[\"NotAvailable\"] when none is found.";

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
Certify::optvalue = "the option name `1` is a symbol with the value `2`, which would replace the name; Clear[`1`] or write the name as the string \"`3`\"";
Certify::optx = "unknown option `1`; the options of Certify are `2`";
LeanVerify::optx = "unknown option `1`; the options of LeanVerify are `2`";
Certify::fracpow = "no Lean certificate is available for `1`: the integrand has a fractional power of the variable, and antiderivatives are certified only for integrands in Q(x)";
Certify::vacuous = "no Lean certificate is available for `1`: the hypothesis `2` != 0 of the Tier 1 statement is zero at every test point under Lean's Real.log (Real.log x = log |x|), so the theorem could be vacuous";
Certify::optignored = "`1` is ignored for a `2` certificate (only identity claims read Domain, Assumptions and Method)";

leanType["Int"] = Integers; leanType["Rat"] = Rationals; leanType["Real"] = Reals;

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
   cannot prove it.  (The Lean encoding: LeanBackend.m, genericStmt.) *)
$genericDomain["CharZero"] = Complexes; $genericDomain["Ordered"] = Reals;
$domainForm["CharZero"] = "CharZeroField"; $domainForm["Ordered"] = "OrderedField";

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
  Module[{vars, st, ty, held, assum, gk, reals, ordered, notReal, dom, sp, auto},
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
  Join[neutralGen[kind, input, result, held,
    <|"Type" -> "Identity", "Statement" -> held, "Variables" -> vars,
      "Pairs" -> {{a, b}}, "Assumptions" -> assum, "Generic" -> gk, "Ring" -> ty|>,
    Join[<|"Form" -> st[[3]], "Variables" -> vars|>, st[[4]]], Join[trusted, st[[5]]],
    st[[2]], {a, b}, If[gk =!= None, $genericDomain[gk], leanType[ty]],
    If[gk =!= None, $domainForm[gk], ty]], <|"Method" -> o["Method"]|>]];

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
(* Sum[t, {k, n}] is Sum[t, {k, 1, n}] (the iterator's definition, as in
   Sum itself), so it is read as that form. *)
cwClaim[Equal[Sum[t_, {k_Symbol, n_}], r_]] := Block[{$cwAnyPow = True},
  inert[Equal, inert[Sum, cw[t], inert[List, csym[k], cw[1], cw[n]]], cw[r]]];
cwClaim[Sum[t_, {k_Symbol, n_}]] := Block[{$cwAnyPow = True},
  inert[Sum, cw[t], inert[List, csym[k], cw[1], cw[n]]]];
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
      If[First[bad] === Power && fracPowIntegrandQ[hc],
        Message[Certify::fracpow, HoldForm[c]],
        Message[Certify::eval, HoldForm[c], First[bad]]];
      Return[$Failed, Module]];
    {cwHeld[t], $cwNotes, $cwTrusted}]];

(* True when an Integrate in the held input has an integrand with a
   non-integer power of a base that depends on the variable (Sqrt[x],
   x^(1/3), ...): outside Q(x), where no certificate exists, rather than a
   claim that "applies Power".  Only the integrand is evaluated, and only to
   choose the message. *)
fracPowIntegrandQ[hc_HoldComplete] := AnyTrue[
  Cases[hc, HoldPattern[Integrate[f_, x_Symbol]] :> {f, x}, {1, 3}],
  !FreeQ[Quiet[#[[1]]], Power[b_, e_] /; !IntegerQ[e] && !FreeQ[b, #[[2]]]] &];

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

(* The neutral generator association (the IR): what was computed and what a
   backend must prove.  "Obligation" holds the held proposition ("Statement")
   and the encoder's data, by "Type"; "Theorem" is what the ProofObject shows.
   No prover text: Lean's "Encode" adds "Statement", "Proof", "Tier", ... *)
neutralGen[kind_, input_, result_, held_, obligation_, cert_, trusted_, hyps_, claim_,
    dom_:Rationals, df_:"Rat"] := <|"Kind" -> kind, "Input" -> input, "Result" -> result,
  "Theorem" -> held, "Obligation" -> obligation, "Certificate" -> cert,
  "TrustedSteps" -> trusted, "Hypotheses" -> hyps, "Domain" -> dom, "Claim" -> claim,
  "DomainForm" -> df, "Method" -> Automatic, "Requires" -> {}|>;

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
  If[o["Tier"] === 1, Return[hasDerivAtIR[input, F, f, x], Module]];
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
  If[!allowedQ[Fs, vars], Return[hasDerivAtIR[input, F, f, x], Module]];
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
  hyps = DeleteDuplicates[Cases[Append[dens, b],
    d_ /; !NumberQ[inertValue[d]] :> (inertHeld[iUnequal[d, 0]] /. back)]];
  neutralGen["Antiderivative", input, F, held,
    <|"Type" -> "PolyIdentities", "Statement" -> held, "Variables" -> vars,
      "Pairs" -> {{lhs, rhs}}, "Ring" -> "Rat"|>,
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

(* HasDerivAt F f x over R (Tier 1, Mathlib): statement, hypotheses and
   witness come from the backend's encoder, replacing these defaults. *)
hasDerivAtIR[input_, F_, f_, x_] := Append[neutralGen["Antiderivative", input, F,
    HoldForm[ForAll[x, D[F, x] == f]],
    <|"Type" -> "HasDerivAt", "Statement" -> HoldForm[ForAll[x, D[F, x] == f]], "F" -> F,
      "f" -> f, "Variable" -> x|>,
    <|"Form" -> "HasDerivAt", "Variable" -> x|>, {}, {}, {D[F, x], f}, Reals, "Real"],
  "Requires" -> {"Mathlib"}];

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
(* The claim's work is done by postEncode, after the backend encoded the
   obligation, so a failed encoding (Certify::tier) is the only message. *)
generate[HoldComplete[Equal[Integrate[f_, x_Symbol], F_]], o_] := Module[{g},
  g = antiderivGen[HoldForm[Integrate[f, x] == F], f, x, o, F];
  If[!AssociationQ[g], Return[g, Module]];
  Append[g, "ClaimForm" -> <|"Type" -> "Antiderivative",
    "Message" -> {HoldForm[Integrate[f, x] == F], HoldForm[D[F, x] == f]},
    "Theorem" -> HoldForm[ForAll[x, D[F, x] == f]],
    "TrustedStep" -> "the claim Integrate[f, x] == F is certified as D[F, x] == f (F is an antiderivative \
of f, shown in \"Theorem\"); Integrate[f, x] == F itself holds only up to an additive term \
with derivative 0, which is locally constant (constant on each interval where F and f are \
defined; it may differ across a point where they are not, as ArcTan[1/x] and -ArcTan[x] \
do at x = 0), not checked"|>]];

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
  neutralGen["Bezout", input, res, held,
    <|"Type" -> "PolyIdentities", "Statement" -> held, "Variables" -> {x},
      "Pairs" -> {{iPlus[iTimes[s, av], iTimes[t, bv]], g}, {av, iTimes[g, qa]},
        {bv, iTimes[g, qb]}}, "Ring" -> ty|>,
    <|"GCD" -> g, "Cofactors" -> {s, t}, "Quotients" -> {qa, qb},
      "PreCheck" -> (Expand[s av + t bv - g] === 0 && Expand[av - g qa] === 0 &&
        Expand[bv - g qb] === 0)|>,
    {}, {}, {{s av + t bv, av, bv}, {g, g qa, g qb}}, leanType[ty], ty]];

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
  neutralGen["IdealMembership", input, res, held,
    <|"Type" -> "PolyIdentities", "Statement" -> held, "Variables" -> vars,
      "Pairs" -> {{qv, iSum[Sequence @@ MapThread[iTimes, {cs, psv}], r]}}, "Ring" -> ty|>,
    <|"Cofactors" -> cs, "Remainder" -> r, "Divisors" -> psv,
      "PreCheck" -> (Expand[qv - cs . psv - r] === 0)|>,
    {"that the remainder is fully reduced with respect to the divisors (not certified; \
the division identity is)"},
    {}, {qv, cs . psv + r}, leanType[ty], ty]];

(* Hypotheses l_i == r_i, goal l == r: q = l - r lies in the ideal of the
   p_i = l_i - r_i.  Proof shape (plan section 0.1): have k_i : p_i = 0,
   have hc : q = Sum c_i * p_i, rw [k_i] at hc, grind. *)
membershipGen[input_, hs_List, {l_, r_}, o_] := Module[{ps, q, vars, cs, eq, held},
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
  neutralGen["IdealMembership", input, cs, held,
    <|"Type" -> "IdealImplication", "Statement" -> held, "Variables" -> vars,
      "Hypotheses" -> hs, "Goal" -> {l, r}, "Cofactors" -> cs|>,
    <|"Cofactors" -> cs, "Generators" -> ps, "Goal" -> q|>, {}, {},
    {q, cs . ps}, $genericDomain["CharZero"], "CharZeroField"]];

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
    held, l, r, alts, sound, dens, eqs, complete = None, soundIR},
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
        complete = neutralGen["Roots", input, sol, held,
          <|"Type" -> "RootsComplete", "Statement" -> held, "Variable" -> x,
            "Equation" -> {l, r}, "Factorization" -> fact|>,
          <|"Completeness" -> True, "Factorization" -> inertValue[fact],
            "Roots" -> roots, "Multiplicities" -> mult|>,
          {}, {}, {Flatten[(l - r) /. sol], 0}]]]];
  (* Soundness: substitute every solution into every constraint (held, so
     the substitution reaches Lean) and let the kernel evaluate.  Lean's
     x / 0 = 0 would accept a "solution" that zeroes a denominator, so each
     denominator of the system is also stated nonzero at each solution.  It
     is also the fallback of the completeness obligation, for a backend that
     cannot state the latter. *)
  dens = denominatorBases[cons[[All, 2 ;; 3]]];
  sound = With[{e = iAnd @@ Join[iUnequal[#, 0] & /@ dens,
       (#[[1]] /. $solveRels) @@ #[[2 ;; 3]] & /@ cons]},
    HoldForm[e]] /. iAnd[c_] :> c;
  held = inertHeld[If[Length[sol] === 1, sound /. First[sol], iAnd @@ ((sound /. #) & /@ sol)]];
  soundIR = neutralGen["Roots", input, sol, held,
    <|"Type" -> "ClosedFact", "Statement" -> held|>,
    <|"Completeness" -> False, "Solutions" -> sol, "Denominators" -> dens|>,
    {"completeness: that Solve returned every solution (only soundness -- each solution \
satisfies each constraint, with every denominator nonzero -- is certified)"},
    {}, {Flatten[Table[(#[[1]] - #[[2]]) & /@ eqs /. s, {s, sol}]], 0}, Rationals, "Closed"];
  If[complete === None, soundIR, Append[complete, "Fallback" -> soundIR]]];

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

(* A closed rational fact about the matrices: the kernel evaluates it. *)
matrixGen[kind_, input_, result_, held_, cert_, trusted_, claim_] :=
  neutralGen[kind, input, result, held, <|"Type" -> "ClosedFact", "Statement" -> held|>,
    cert, trusted, {}, claim, Rationals, "Closed"];

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
    st, terms, held, used},
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
  held = inertHeld[If[conds === {}, iForAll[vars, inertRel[rel, a, b]],
    iForAll[vars, If[Length[conds] === 1, First[conds], iAnd @@ conds], inertRel[rel, a, b]]]];
  neutralGen["Nonnegativity", input, True, held,
    <|"Type" -> "WeightedSOS", "Statement" -> held, "Variables" -> vars,
      "Conditions" -> conds, "Relation" -> rel, "Target" -> target, "Squares" -> terms,
      "Multipliers" -> used|>,
    <|"Squares" -> terms, "Method" -> st[[2]], "Target" -> target, "Multipliers" -> used,
      "PreCheck" -> True|>,
    {}, {}, None, $genericDomain["Ordered"], "OrderedField"]];

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
    sigma, st, terms, held, rel, feas, value, lamFull},
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
  neutralGen["Optimality", input, res, held,
    <|"Type" -> "KKT", "Statement" -> held, "Variables" -> vars, "Objective" -> f,
      "Value" -> m, "Minimizer" -> xs, "Sign" -> sgn, "Constraints" -> cons,
      "Eliminated" -> eqHyps, "Reduced" -> gr, "Multipliers" -> lamFull, "Squares" -> terms,
      "Gap" -> Expand[Fr - M], "Facts" -> Length[feas] + 1|>,
    <|"Value" -> m, "Minimizer" -> xs, "Multipliers" -> lamFull,
      "Constraints" -> ineqs[[All, 3]] /. $inertRules, "Eliminated" -> (Rule @@@ eqHyps),
      "Squares" -> terms, "Method" -> st[[2]], "LinearProgram" -> (sigma === 0),
      "PreCheck" -> (Expand[F /. xs] === M && AllTrue[ineqs[[All, 2]] /. xs, # >= 0 &] &&
        AllTrue[Select[cons, #[[1]] === "eq" &][[All, 2]] /. xs, # === 0 &])|>,
    {}, {}, None, $genericDomain["Ordered"], "OrderedField"]];

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

generate[HoldComplete[FactorInteger[n_]], o_] := Module[{input = HoldForm[FactorInteger[n]],
    nv = n, fi, odd, certs, held},
  If[!IntegerQ[nv] || nv < 2, Return[nocert[input], Module]];
  fi = FactorInteger[nv];
  odd = DeleteCases[fi[[All, 1]], 2];
  certs = pratt /@ odd;
  If[MemberQ[certs, $Failed], Return[nocert[input], Module]];
  held = inertHeld[iAnd @@ Join[{iEqual[nv, iProd @@ (If[#[[2]] === 1, #[[1]], iPower @@ #] & /@ fi)]},
    iPrimeQ /@ fi[[All, 1]]]];
  neutralGen["Primality", input, fi, held, primalityOb["Factorization", nv, held,
      "Factorization" -> fi, "Certificates" -> certs],
    <|"Factorization" -> fi, "Pratt" -> AssociationThread[odd -> certs]|>,
    {$lucasStep},
    {}, {Times @@ (#[[1]]^#[[2]] & /@ fi), nv}, Integers, "Closed"]];

(* The "Primality" obligation: n is prime (a flat Pratt certificate per odd
   prime), 2 is prime, n < 2 or a composite n is not prime. *)
primalityOb[form_, nv_, held_, data___Rule] := <|"Type" -> "Primality", "Statement" -> held,
  "Form" -> form, "N" -> nv, data|>;

generate[HoldComplete[PrimeQ[n_]], o_] := Module[{input = HoldForm[PrimeQ[n]], nv = n, a, cert},
  If[!IntegerQ[nv] || nv < 0, Return[nocert[input], Module]];
  Which[
    nv < 2,
      With[{h = With[{x = nv}, HoldForm[!PrimeQ[x]]]},
        neutralGen["Primality", input, False, h, primalityOb["Small", nv, h],
          <|"Witness" -> "n < 2"|>, {}, {}, None, Integers, "Closed"]],
    nv === 2,
      neutralGen["Primality", input, True, HoldForm[PrimeQ[2]],
        primalityOb["Two", nv, HoldForm[PrimeQ[2]]],
        <|"Witness" -> "the divisors of 2"|>, {}, {}, None, Integers, "Closed"],
    PrimeQ[nv],
      cert = pratt[nv];
      If[cert === $Failed, Return[nocert[input], Module]];
      With[{h = With[{x = nv}, HoldForm[PrimeQ[x]]]},
        neutralGen["Primality", input, True, h,
          primalityOb["Prime", nv, h, "Certificates" -> {cert}],
          <|"Pratt" -> cert|>, {$lucasStep}, {}, None, Integers, "Closed"]],
    True,
      a = FactorInteger[nv][[1, 1]];
      With[{h = With[{x = nv}, HoldForm[!PrimeQ[x]]]},
        neutralGen["Primality", input, False, h,
          primalityOb["Composite", nv, h, "Divisor" -> a],
          <|"Factor" -> {a, nv/a}|>, {}, {}, None, Integers, "Closed"]]]];

(* ---------------------------------------------------------------------- *)
(* Row 16: Sum[t, {k, a, n}] with t = polynomial(k) c^k (c rational) or a   *)
(* polynomial: Lean defines t and S by recursion from a and proves          *)
(* forall m : Nat, S m = closed(m + a) by induction.                        *)
(* ---------------------------------------------------------------------- *)

sumGen[input_, t_, k_Symbol, a_, n_Symbol, cl0_, o_] := Module[{av, held, closed},
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
  (* Inert heads, not With: With renames the iterator of a held Sum. *)
  held = inertHeld[iForAll[n, iGreaterEqual[n, av], iEqual[iSumOp[t, {k, av, n}], closed]]];
  neutralGen["ClosedFormSum", input, closed, held,
    <|"Type" -> "Induction", "Statement" -> held, "Term" -> t, "Variable" -> k,
      "LowerLimit" -> av, "SumVariable" -> n, "ClosedForm" -> closed|>,
    <|"Term" -> t, "Variable" -> k, "LowerLimit" -> av, "ClosedForm" -> closed,
      "PreCheck" -> (Expand[(closed /. n -> av) - (t /. k -> av)] === 0 &&
        Expand[(closed /. n -> n + 1) - closed - (t /. k -> n + 1)] === 0)|>,
    {"the closed form is proved for integer n >= " <> ToString[av] <> " only (Lean's \
statement is over n = " <> If[av === 0, "", ToString[av] <> " + "] <> "m with m : Nat), as \
\"Theorem\" states"},
    {}, {Sum @@ {t, {k, av, n}}, closed}, Rationals, "Int"]];

(* The held input is passed whole (s, c): wrapping a held Sum in With would
   rename its iterator. *)
generate[HoldComplete[s:Sum[t_, {k_Symbol, a_, n_Symbol}]], o_] :=
  sumGen[HoldForm[s], t, k, a, n, Automatic, o];

generate[HoldComplete[c:Equal[Sum[t_, {k_Symbol, a_, n_Symbol}], rhs_]], o_] :=
  sumGen[HoldForm[c], t, k, a, n, rhs, o];

(* Anything else: a symbol or call holding a supported shape (c = a == b;
   Certify[c]) is evaluated once and dispatched again; a fixed point that
   matches no rule has no certificate. *)
(* A supported computation head in a form no generator above takes (an
   extra argument, as in Solve[eqs, x, Reals]) has no certificate; it is
   not evaluated and re-dispatched, so the message names the input as
   written, not its value. *)
generate[HoldComplete[c:(h_Symbol)[___]], o_] /; MemberQ[$cwComputations, h] :=
  nocert[HoldForm[c]];

generate[HoldComplete[c_], o_] := With[{v = c},
  If[HoldComplete[v] =!= HoldComplete[c], generate[HoldComplete[v], o],
    Message[Certify::nocert, HoldForm[c]]; $Failed]];

(* ---------------------------------------------------------------------- *)
(* Certify / ProofObject                                                   *)
(* ---------------------------------------------------------------------- *)

Options[Certify] = {"Engine" -> "Lean", "Tier" -> Automatic, "Domain" -> Automatic,
  Method -> Automatic, TimeConstraint -> 60, Assumptions -> {}, "Verify" -> True,
  "LeanExecutable" -> Automatic, "LeanProject" -> Automatic};

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
  (* An option name that is a symbol with a value would be replaced by that
     value when the options are evaluated (Engine = "zzz" turns Engine -> v
     into "zzz" -> v): refuse it while the name is still held. *)
  bad = Cases[HoldComplete[opts], (Rule | RuleDelayed)[s_Symbol, _] /;
      OwnValues[s] =!= {} :> {HoldForm[s], ToString[OwnValues[s][[1, 2]], InputForm], SymbolName[Unevaluated[s]]},
    {1, 2}];
  If[bad =!= {},
    Message[Certify::optvalue, Sequence @@ bad[[1]]];
    Return[$Failed, Module]];
  as = Cases[HoldComplete[opts], (Rule | RuleDelayed)[Assumptions | "Assumptions", a_] :>
    HoldComplete[a], {1, 2}];
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

certOptions[opts_List] := Module[{keys = {"Engine", "Tier", "Domain", "Method",
    "TimeConstraint", "Assumptions", "Verify", "LeanExecutable", "LeanProject"}},
  If[optionNameProblem[Certify, opts] =!= None, Return[$Failed, Module]];
  AssociationThread[keys, OptionValue[Certify, opts, #] & /@ keys]];

(* None, or $Failed after f::optx, when an evaluated option's name is not
   one of Options[f] (compared by text, so a symbol or a string of that name
   is accepted).  A typo such as Teir -> 1 must not be ignored silently. *)
optionNameProblem[f_, opts_List] := Module[{known = optText /@ Keys[Options[f]], bad},
  bad = Select[Flatten[opts], !MatchQ[#, (Rule | RuleDelayed)[n_, _] /;
      MemberQ[known, optText[n]]] &];
  If[bad === {}, None,
    Message[MessageName[f, "optx"], bad[[1]], StringRiffle[known, ", "]]; $Failed]];
optText[s_String] := s;
optText[s_Symbol] := SymbolName[Unevaluated[s]];
optText[_] := None;

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

(* ---------------------------------------------------------------------- *)
(* Backends and the pipeline                                               *)
(* ---------------------------------------------------------------------- *)

(* The backend registry: name -> spec, registered at load by the file that
   defines the backend.  Spec keys: "Name", "Trust" (Kernel | ExactArithmetic
   | Heuristic), "Encode" ([ir, o] -> artefacts | $Failed after a message |
   "Untranslatable" when the IR has a "Fallback"), "Run" ([items, o] ->
   verdicts; item <|"Index", "G"|>), "Present" ([item, v, o] -> ProofObject),
   "Label" ([assoc] -> the summary box's third field). *)
$backends = <||>;
$backendKeys = {"Name", "Trust", "Encode", "Run", "Present", "Label"};
registerBackend[name_String, spec_Association] := If[
  SubsetQ[Keys[spec], $backendKeys] && spec["Name"] === name &&
    MemberQ[{"Kernel", "ExactArithmetic", "Heuristic"}, spec["Trust"]],
  $backends = Append[$backends, name -> spec]; name, $Failed];

(* The IR joined with the engine's artefacts, or $Failed after a message; run
   right after the generator, so an encoding refusal keeps its position. *)
encodeFor[ir_, engine_, o_] := With[{art = $backends[engine]["Encode"][ir, o]},
  Which[
    art === "Untranslatable" && KeyExistsQ[ir, "Fallback"], encodeFor[ir["Fallback"], engine, o],
    art === "Untranslatable", nocert[ir["Input"]],
    AssociationQ[art], Join[ir, art],
    True, $Failed]];

(* The claim-form work deferred until the obligation is encoded: the claim
   Integrate[f, x] == F is certified as D[F, x] == f (Certify::antideriv,
   the Tier 0 "Theorem", the trusted step on the gap). *)
postEncode[g_] := If[!KeyExistsQ[g, "ClaimForm"], g, With[{cf = g["ClaimForm"]},
  Message[Certify::antideriv, Sequence @@ cf["Message"]];
  Join[g, <|"Theorem" -> If[g["Tier"] === 1, g["Theorem"], cf["Theorem"]],
    "TrustedSteps" -> Append[g["TrustedSteps"], cf["TrustedStep"]]|>]]];

(* The generator association of one item plus the tier it runs at
   ("TierUsed"), or $Failed (after a message). *)
prepare[hc_HoldComplete, o_] := Module[{g, tier, proj, rw, sides, ign},
  (* The input as written (rewriteInput); Certify::eval when it applies a
     computation Mathilda would have to evaluate first. *)
  rw = rewriteInput[hc];
  If[rw === $Failed, Return[$Failed, Module]];
  g = generate[First[rw], o];
  If[!AssociationQ[g], Return[$Failed, Module]];
  g = encodeFor[g, o["Engine"], o];
  If[!AssociationQ[g], Return[$Failed, Module]];
  g = postEncode[g];
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
    ign = Select[{"Domain" -> o["Domain"], Assumptions -> o["Assumptions"],
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

(* The heuristic shows the Lean statement and batch source. *)
heuristicRun[items_List, o_] := MapThread[Append[heuristicVerdict[#1["G"]], "Source" -> #2] &,
  {items, leanBatchSources[items, o]}];
heuristicPresent[item_, v_, o_] := With[{g = leanNamed[item]},
  mkProofObject[Append[g, "TrustedSteps" ->
      {"Mathilda's PossibleZeroQ / Simplify (a heuristic check, not a proof)"}],
    v, "Mathilda", "Heuristic", Missing["NotChecked"], v["Source"], g["Statement"]]];

registerBackend["Mathilda", <|"Name" -> "Mathilda", "Trust" -> "Heuristic",
  "Encode" -> leanEncode, "Run" -> heuristicRun, "Present" -> heuristicPresent,
  "Label" -> ("Mathilda Heuristic" &), "Options" -> {}|>];

(* Certify a list of held items with the engine's backend (Lean: one run per
   tier), aligned with the input ($Failed: no certificate).  Only a kernel
   backend's True is shown as True. *)
certifyItems[hcs_List, opts_List] := Module[{o, engine, gens, items, vs, v, out},
  o = certOptions[opts];
  If[o === $Failed, Return[$Failed, Module]];
  engine = o["Engine"];
  If[!StringQ[engine] || !KeyExistsQ[$backends, engine],
    Message[Certify::engine, engine]; Return[$Failed, Module]];
  (* Every Domain value has one stated meaning ($domains); anything else
     (Algebraics, "Complexes", ...) is refused rather than silently mapped. *)
  If[!AnyTrue[$domains, # === o["Domain"] &],
    Message[Certify::domain, ToString[o["Domain"], InputForm]]; Return[$Failed, Module]];
  With[{why = methodProblem[o["Method"]]},
    If[why =!= None, Message[Certify::method, o["Method"], why]; Return[$Failed, Module]]];
  gens = prepare[#, o] & /@ hcs;
  items = Table[<|"Index" -> k, "G" -> gens[[k]]|>,
    {k, Select[Range[Length[gens]], AssociationQ[gens[[#]]] &]}];
  out = ConstantArray[$Failed, Length[gens]];
  If[items =!= {},
    vs = $backends[engine]["Run"][items, o];
    Do[
      v = vs[[j]];
      If[v["Verified"] === True && $backends[engine]["Trust"] =!= "Kernel",
        v = Append[v, "Verified" -> Missing[engine, True]]];
      out[[items[[j]]["Index"]]] = $backends[engine]["Present"][items[[j]], v, o],
      {j, Length[items]}]];
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
  With[{e = a["Engine"]}, If[StringQ[e] && KeyExistsQ[$backends, e],
    $backends[e]["Label"][a], "Mathilda Heuristic"]],
  " | ", timingText[a["Timing"]], "]"];

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

End[];

EndPackage[];

(* The backends (each opens Lean`Private` itself; no Lean` symbol). *)
LoadModule["lean/LeanBackend.m"];

SetAttributes[{Lean`Certify, Lean`ProofObject, Lean`LeanVerify}, Protected];
