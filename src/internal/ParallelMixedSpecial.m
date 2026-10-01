(* ParallelMixedSpecial.m -- see the header comment below.
   -------------------------------------------------------------------------
   MATHILDA LOADING CONTRACT.  This file is NOT a package: it carries no
   BeginPackage and no Begin["`Private`"], and it is loaded INTO
   ParallelMixed`Private` -- the same mechanism mixed/logrewrite.m uses.

   The research original is a package (BeginPackage["ParallelMixedSpecial`",
   {"ParallelMixed`"}]) that then does

       AppendTo[$ContextPath, "ParallelMixed`Private`"];

   so that the 41 internals of Part II it needs (Can, TMul, TDiv, PlaceBounds,
   SamplePoint, DecideBound, $analyses, ...) resolve by their short names.  That
   does not work here: Mathilda resolves a short private *variable* through
   $ContextPath but not a short private *function*, so every one of those 41
   names would have stayed unevaluated.  Being loaded into that context instead
   makes the short names resolve natively, because they really are in it -- and
   it costs nothing, since this stage is an extension of Part II rather than an
   independent package.

   Consequences, all deliberate:
     * the eight public entry points are ParallelMixed`Private` symbols.  The C
       method addresses IntegrateSurfaceSpecial by its qualified name (see
       builtin_integrate_pms in src/calculus/integrate.c), exactly as the
       ParallelMixedTower method addresses ParallelMixed`ParallelIntegrateMixed.
     * the ::usage strings below are kept: `f::usage = "..."` also registers the
       string as f's docstring, so ?name works on them.
     * Part II is loaded through LoadModule rather than the original's
       Get[FileNameJoin[{DirectoryName[$InputFileName], ...}]] -- $InputFileName
       and DirectoryName are not implemented (MATHILDA_DIVERGENCES.md A22), and
       LoadModule additionally carries the load-once memo that stops the 217 KB
       Part II file being evaluated twice in one session.

   Divergences from the research original, each forced and each noted at its
   site: FreeQ against a bare type symbol does not match an atom's type-head, so
   the five tests for a Gaussian constant use the pattern _Complex (A19).  The
   iterator-capture renames that mixed/ParallelMixed.m needs are NOT needed here:
   A18 was fixed at the root in v0.238, so Do, Sum, Product and Table are all
   capture-avoiding and this file keeps the original's iterator names.
   ------------------------------------------------------------------------- *)

(* Part II.  Idempotent: LoadModule is memoised on the relative path. *)
If[! TrueQ[LoadModule["mixed/ParallelMixed.m"]],
  Print["ParallelMixedSpecial: mixed/ParallelMixed.m not found under src/internal; \
the special-function stage cannot run"]];

(* Everything below is evaluated INSIDE ParallelMixed`Private` (see above), so
   Part II's internals are in scope by their short names. *)
Begin["ParallelMixed`Private`"];

(* ParallelMixedSpecial.m -- the Mathilda port of parallel_mixed_special.py

   Parallel integration over mixed towers with antiderivatives in special
   functions (paper: rn-radicals-specfun.tex, boxes S1--S10), on top of the
   Part II port ParallelMixed.m (mixed/ParallelMixed.m, loaded below).

   The class: elementary functions and the upper incomplete gamma function
   Gamma[s, .] at rational s in [0, 1) -- ExpIntegralEi, LogIntegral,
   SinIntegral, CosIntegral (s = 0), Erf, Erfc, Erfi (s = 1/2), Gamma[1/m, .]
   in general -- and the incomplete elliptic integrals EllipticF, EllipticE,
   EllipticPi of constant modulus carried by a radical y^2 = q(g) whose radicand
   is a polynomial in one generator g of the tower (a pencil).  Every special
   function enters through a KERNEL -- a known element of the tower with a
   known antiderivative -- as one more column of Part II's linear system with
   an unknown constant coefficient, or, when the residues determine its
   coefficient, as a known subtraction before the system is assembled.

   Box to function map (the Python names in CamelCase):
     S1  ExpSources, InferSources          S6  SolveWithKernels, Minimal
     S2  EiFromResidues                    S7  ParallelIntegrateSpecial (Pis)
     S3  GammaCandidates                   S8  CertifyNonelementary, IndependenceCertificate
     S4  EllipticColumns, ThirdKind        S9  Present, KernelAntiderivative, FixBranch
     S5  DecideBoundExt, SpecialDataExt    S10 PartialIntegrate
         (ExtendedBounds)
   Entry points: IntegrateSurfaceSpecial[f, x], IntegrateSurfacePartial[f, x],
   ParallelIntegrateSpecial[{f0, f1}, T], PartialIntegrate[{f0, f1}, T].

   The Rothstein--Caviness structure step of special/build_tower.py is
   BuildTower[..., "StructureTheorem" -> True] in ParallelMixed.wl; the entry
   points here turn it on.  The bound extensions (K5)--(K8), (B) are installed
   by ExtendedBounds, a Block over the four hooks DecideBound, SpecialData,
   InfConsts, SpecialConsts of ParallelMixed` (the originals stay callable as
   DecideBoundOrig, ...), the counterpart of the monkeypatching of
   extended_bounds.

   Representations: a kernel is an Association <|"kind" -> "gamma" | "ell" |
   "third", "col" -> tuple, ...|>; a special answer <|"type" -> "special",
   "elementary", "terms" -> {{c, K}, ...}, "certified", "certificate",
   "independence", "necessity", "base"|>; a partial answer <|"type" ->
   "partial", "elementary", "terms", "remainder", "remainderCertified"|>; a
   status is a list {"failed", ...}, {"not in class", ...}, as in Part II; an
   elementary answer is an expression.  The unevaluated integral is
   Inactive[Integrate].  The verbose trace lines are those of the Python
   module.                                                                   *)

IntegrateSurfaceSpecial::usage = "IntegrateSurfaceSpecial[f, x, opts] integrates the expression f in x in the class of elementary functions, Gamma[s, .] at rational s (ExpIntegralEi, LogIntegral, Erf, Erfc, Erfi, SinIntegral, CosIntegral, Gamma[1/m, .]) and elliptic integrals of a pencil. The tower is built by BuildTower with the structure theorem (or given by \"Tower\" -> {T, {f0, f1}, back}); returns {answer, verified} or a status list {\"failed\", ...} / {\"not elementary\", ...} / {\"not in class\", ...}. Options: \"Verbose\", \"Tower\", \"Samples\", \"Strict\" (default True: a special answer only with a certificate of non-elementarity), \"Details\".";
IntegrateSurfacePartial::usage = "IntegrateSurfacePartial[f, x, opts] is IntegrateSurfaceSpecial with a partial answer I + Inactive[Integrate][r, x] when no complete one is found (\"Special\" -> False: the elementary partial answer, no special function at all). Returns {answer, verified} or a status.";
ParallelIntegrateSpecial::usage = "ParallelIntegrateSpecial[{f0, f1}, T, opts] integrates f0 + f1 y over the tower T with special functions (Algorithm S7): an expression (elementary), a special result Association, or a status list.";
PartialIntegrate::usage = "PartialIntegrate[{f0, f1}, T, opts] is Algorithm S10: a complete answer if ParallelIntegrateSpecial finds one, else a partial result Association f = D(I) + r with a reduced remainder r.";
ExpSources::usage = "ExpSources[T, back, Y, verbose] lists the exponential sources {name, xi, lambda} of the tower (Algorithm S1).";
Present::usage = "Present[e] rewrites Gamma[1/2, .] through Erfc, Erfc[I z] through Erfi, ExpIntegralEi[Log[a]] through LogIntegral and conjugate ExpIntegralEi pairs through CosIntegral and SinIntegral.";
KernelAntiderivative::usage = "KernelAntiderivative[K] is the antiderivative of the kernel K in the tower generators.";
ExtendedBounds::usage = "ExtendedBounds[on, body] evaluates body with the bound criteria (K5)--(K8) and the branch places (B) of Section 5 installed into Part II's Algorithm 6.";

$ExtendedBoundsSP = True;       (* (K5)-(K8) and (B) of Section 5 in the special stage *)
$StrictSP = True;               (* Theorem 8.6: special answers only with a certificate of non-elementarity *)
$spCounter = 0;

(* ----------------------------------------------------------- stage profiling *)
(* $SpecialProfile -> True makes each pipeline stage below accumulate its wall
   time and its call count into $SpecialStages; SpecialProfile[] reads it back,
   sorted by cost.  This exists because a `sample` profile cannot attribute this
   pipeline at all: every .m stage bottoms out in the same generic C primitives
   (evaluate_step, the matcher, malloc churn), so the profile names those and
   not the stage.  The stages are millisecond-scale, so the TrueQ guard on the
   off path costs nothing measurable.  AbsoluteTiming, not AbsoluteTime:
   AbsoluteTime[] has integer-second resolution here. *)
$SpecialProfile = False;
$SpecialStages = <||>;
SetAttributes[PMStage, HoldRest];
PMStage[lbl_, body_] := If[! TrueQ[$SpecialProfile], body,
  Module[{r = AbsoluteTiming[body], p},
    p = Lookup[$SpecialStages, lbl, {0., 0}];
    $SpecialStages[lbl] = {p[[1]] + r[[1]], p[[2]] + 1};
    r[[2]]]];
SpecialProfile[] := Reverse[SortBy[Normal[$SpecialStages], #[[2, 1]] &]];
SpecialProfileReset[] := ($SpecialStages = <||>;);

(* Robustness backstop, the counterpart of Part II's $ParallelMixedTimeBudget
   (mixed/ParallelMixed.m) and needed for the same reason: the kernel search of
   Algorithm S6 escalates through six {split, retry} configurations, and on an
   integrand with two transcendental generators and no special kernel that can
   close it the ansatz grows until the linear solve dominates -- unboundedly, as
   interactive Integrate has no timeout of its own.  Without this, one such
   integrand turns a plain Integrate[] into a minutes-long grind
   (Sin[x^2 + Log[x]] Cos[x] measured >120 s).  A wall-clock budget converts it
   into a clean {"failed", ...} decline, which the cascade and the Method
   surface both already handle.  Set for parity with Part II; deterministic
   ansatz-size caps are follow-up work, as they are there.

   The Automatic cascade scopes this DOWN to PMS_CASCADE_BUDGET_SECONDS (10 s;
   src/calculus/integrate.c, which Blocks this symbol -- UNQUALIFIED, since the
   context injection leaves it in Global`) : a stage applied to every integrand
   nothing else closed must be cheap, where an explicit
   Method -> "ParallelMixedSpecial" is a deliberate request for the full search. *)
$SpecialTimeBudget = 45;

(* ------------------------------------------------------------------ helpers *)
ToStr[e_] := ToString[e, InputForm];
FreeSyms[e_] := DeleteDuplicates[Cases[{e}, s_Symbol /; ! NumericQ[s] &&
    ! MemberQ[{ComplexInfinity, Indeterminate, None, True, False, Null, Infinity}, s], {0, Infinity}]];
HasAny[e_, syms_List] := syms =!= {} && ! FreeQ[e, Alternatives @@ syms];
NCo[T_] := If[T["q"] === None, 1, T["n"]];
TZero[T_] := ConstantArray[0, T["n"]];
AllZeroQ[u_] := AllTrue[u, Can[#] === 0 &];
FiniteValQ[v_] := FreeQ[v, ComplexInfinity | Indeterminate | DirectedInfinity | Infinity];
CanT[u_] := Can /@ u;
(* CanR[e]: Can, with a denominator that Cancel over the extension left factored over
   Q(i) although the quotient has rational coefficients multiplied out again (the
   form of Python's _c: t^2 + 1, not (t - I)(t + I)), so that Part II's analysis
   factors it over Q *)
CanR[e_] := With[{c = Can[e]}, If[FreeQ[c, _Complex], c,
  With[{r = Together[Expand[Numerator[c]]/Expand[Denominator[c]]]}, If[FreeQ[r, _Complex], Can[r], c]]]];
CanRT[u_] := CanR /@ u;
SetAttributes[QuietCheck, HoldAll];
QuietCheck[e_, fail_] := Quiet[Check[e, fail]];
(* exact zero test of a constant or an expression *)
ExactZeroQ[z_] := With[{w = Quiet[Together[z]]}, Which[w === 0, True,
  NumericQ[w], TrueQ[Quiet[RootReduce[w] === 0]] || Quiet[Simplify[w]] === 0,
  True, Quiet[Simplify[w]] === 0]];
(* the irreducible factors of a polynomial (sympy's factor_list(Poly(den, *gens))[1]),
   over the Gaussian rationals when the polynomial has the constant I *)
FactorsOf[poly_] := Module[{fl},
  fl = If[FreeQ[poly, _Complex], FactorList[poly], Quiet[FactorList[poly, GaussianIntegers -> True]]];
  Select[fl, ! NumericQ[#[[1]]] &]];
(* the roots in radicals of a polynomial in g (sympy's roots: Root objects of an
   unsolvable factor are left out), with multiplicities: {{root, mult}, ...} *)
RootsMult[poly_, g_] := Module[{out = {}, sols},
  Do[If[FreeQ[fac[[1]], g], Continue[]];
    sols = Quiet[g /. Solve[fac[[1]] == 0, g, Cubics -> True, Quartics -> True]];
    If[! ListQ[sols], Continue[]];
    Do[If[FreeQ[r, Root], AppendTo[out, {r, fac[[2]]}]], {r, DeleteDuplicates[sols]}],
    {fac, FactorList[poly]}];
  out];
RootsOf[poly_, g_] := RootsMult[poly, g][[All, 1]];
(* numerator coefficients of an expression as a polynomial in the generators *)
NumCoeffs[e_, gens_] := Module[{num = Expand[Numerator[Together[e]]], rest},
  rest = Select[gens, ! FreeQ[num, #] &];
  Which[num === 0, {}, rest === {}, {num}, True, CoefficientRules[num, rest][[All, 2]]]];
(* LinSolveZero[eqs, vars]: sympy's linsolve with the free unknowns set to 0 -- the
   row-reduced (column order vars) particular solution; None if inconsistent *)
LinSolveZero[eqs_, vars_] := Module[{b, M, aug, red, n = Length[vars], xs, ok = True, p},
  If[eqs === {}, Return[ConstantArray[0, n]]];
  {b, M} = Normal[CoefficientArrays[eqs, vars]];
  aug = MapThread[Append, {M, -b}];
  red = RowReduce[aug, ZeroTest -> ExactZeroQ];
  xs = ConstantArray[0, n];
  Do[p = LengthWhile[red[[k]], ExactZeroQ] + 1;
    Which[p > n + 1, Null, p == n + 1, ok = False, True, xs[[p]] = red[[k, n + 1]]],
    {k, Length[red]}];
  If[ok, xs, None]];
(* sympy's could_extract_minus_sign, syntactically *)
MinusSignQ[z_] := Internal`SyntacticNegativeQ[z] || (Head[z] === Plus && Internal`SyntacticNegativeQ[First[z]]);
(* count_ops in the manner of sympy (the operations of the expression tree) *)
CountOps[e_] := Which[
  NumericQ[e] && AtomQ[e], Which[Head[e] === Rational, 1 + If[e < 0, 1, 0], Head[e] === Integer, If[e < 0, 1, 0], True, 0],
  AtomQ[e], 0,
  Head[e] === Plus, Length[e] - 1 + Total[CountOps /@ (List @@ e)],
  Head[e] === Times, Module[{num = Numerator[e], den = Denominator[e], fn, fd},
    fn = If[Head[num] === Times, List @@ num, {num}]; fd = If[Head[den] === Times, List @@ den, {den}];
    fn = DeleteCases[fn, 1]; fd = DeleteCases[fd, 1];
    Max[0, Length[fn] - 1] + If[fd === {}, 0, Length[fd]] + Total[CountOps /@ Join[Abs /@ fn, fd]] +
      If[Internal`SyntacticNegativeQ[e], 1, 0]],
  Head[e] === Power, If[e[[2]] === -1, 1 + CountOps[e[[1]]], 1 + CountOps[e[[1]]] + CountOps[e[[2]]]],
  True, 1 + Total[CountOps /@ (List @@ e)]];

(* -------------------------------------------- additive decomposition (S11) *)
(* The cost of the stage is MULTIPLICATIVE in the number of independent
   generator families the integrand mentions, because one ansatz is built
   spanning the product of the families (AnsatzSystem, ParallelMixed.m:2966,
   driven by the {split, retry} ladder below and re-run per drop by Minimal).
   Measured on this binary, v0.244:

     Log[x]/x                        0.021 s
     Sin[x]/x                        0.216 s
     Exp[-x^2]                       0.063 s
     Log[x]/x + Sin[x]/x             5.75 s   (24x the sum of the parts)
     Log[x]/x + Sin[x]/x + Exp[-x^2] 20.12 s  (67x)

   and in each case the answer is the concatenation of the independent
   answers.  So: split the integrand into additive blocks whose generator
   families are pairwise disjoint, integrate each block on its own tower, and
   sum.  Per-block BuildTower is ~4-11 ms, noise against the per-call floor,
   and it is strictly MORE capable than the joint one -- two different
   radicals in two different terms (Sqrt[1-x^2] + Sqrt[1+x^3]) make the joint
   tower fail outright and each block succeed.

   SOUNDNESS.  Disjoint families do NOT by themselves make the split
   complete.  The obvious argument -- "terms that cancel must share a
   generator" -- is FALSE, and this counterexample is why the guard below
   exists (both halves verified on this binary):

     Log[1+x]/x + Log[x]/(1+x)  ->  Log[x] Log[1+x]   elementary
     Log[1+x]/x   alone         ->  "not in class"
     Log[x]/(1+x) alone         ->  "not in class"

   The families {Log[1+x]} and {Log[x]} are disjoint, the sum is elementary,
   and neither half is even in the class.  The mechanism is specific to
   PRIMITIVES: if t_A, t_B are primitives whose derivatives lie in the shared
   base C(x), then

       D(c t_A t_B) = c (D t_A) t_B + c t_A (D t_B)

   is a two-term sum whose terms have DISJOINT support, so Liouville's v_0 is
   entitled to the mixed monomial t_A t_B.  No other generator kind can do
   this: a hyperexponential has D theta = lambda theta and so keeps theta in
   every term; likewise Tan/Tanh (D t = (1 +- t^2) Da), ProductLog, and a
   flattened root.  An InvRadical primitive (ArcSin, ArcSec, ...) carries its
   radical into D t, which is not in the shared base -- and if the radical
   WERE shared the two blocks would already have merged on it.

   Writing L = L_B(t_A) and v_0 = Sum_k v_k t_A^k with v_k in L_B, the
   t_A^n coefficient for the top degree n >= 1 gives D v_n in C(x): block B
   must itself hold an element whose derivative is in the shared base.  Call
   such a block COUPLING-CAPABLE (BasePrimitiveQ below: it carries a Log or an
   InvRational head at an argument rational in x alone).  Hence

     Lemma.  If f = Sum_B f_B over blocks with pairwise-disjoint families and
     Int f is elementary while some Int f_B is not, then at least TWO blocks
     are coupling-capable and the obstruction is a constant-coefficient
     bilinear form in base primitives drawn from distinct blocks.

   Two corollaries drive the guard: with at most one coupling-capable block
   among those that did not close, no cross-block rescue exists; and a block
   that closes completely cannot take part in a coupling at all.  The residual
   assumption is the CLASS boundary rather than the elementary one -- a kernel
   built from a joint exponential monomial (JointMonomials) has mixed support,
   and I have neither a proof that it cannot rescue a block-pure sum nor an
   integrand where it does.  Completeness here is already heuristic (the six
   {split, retry} rungs, the guessed bounds, the budget decline), and answer
   soundness does not rest on it.                                           *)

(* a^b with a non-rational exponent is exp(b log a), held as PMExp so the
   kernel does not fold it back: the same rewrite BuildTower applies at
   ParallelMixed.m:2424, so the scan sees x^Sqrt[2] and 2^x as exponentials *)
PMExpHold[e_, x_] := e /. Power[b_, ex_] /; b =!= E && ! MatchQ[ex, _Rational | _Integer] &&
  (! FreeQ[ex, x] || ! FreeQ[b, x]) :> PMExp[ex Log[b]];

(* the rational content divided out.  Exp[x], Exp[-x] and Exp[2 x] are ONE
   family (one eulerStep generator, and the structure theorem at
   ParallelMixed.m:2466 refines Exp[x/3] and Exp[3x/2] to a common e^(x/6)),
   as are Tan[x/2] and Tan[x]; Exp[x] and Exp[x^2] are two.
   Only RATIONAL content is divided out, not every numeric factor: 2^x is
   PMExp[x Log[2]], and e^(x Log 2) is NOT algebraic over e^x (Log[2] is
   irrational), so stripping Log[2] would wrongly fuse 2^x with Exp[x]. *)
RatContent[a_] := Which[a === 0, 1, MatchQ[a, _Integer | _Rational], a,
  Head[a] === Times,
    With[{c = Times @@ Select[List @@ a, MatchQ[#, _Integer | _Rational] &]},
      If[MatchQ[c, _Integer | _Rational] && c =!= 0, c, 1]],
  True, 1];
NormTerm[a_] := With[{c = RatContent[a]}, If[c =!= 0, Expand[a/c], a]];
(* a polynomial with its integer content divided out, so that the radicands
   of Sqrt[4x^2+4] and Sqrt[x^2+1] are one family.  The SIGN is deliberately
   left alone: Sqrt[x-x^3] and Sqrt[-x+x^3] are different radicands. *)
NormPoly[a_] := Module[{e = Quiet[Expand[Together[a]]], ft},
  ft = Quiet[Check[FactorTermsList[e], {1, e}]];
  NormTerm[If[ListQ[ft] && Length[ft] == 2, ft[[2]], e]]];
(* an exponent contributes one key per additive term, because BuildTower's
   exponential split (ParallelMixed.m:2449) breaks e^(x - t) into e^x e^(-t);
   purely numeric terms are a constant factor and carry no family *)
ExpKeys[a_] := Module[{t = Together[a], ts},
  ts = If[Head[Expand[t]] === Plus, List @@ Expand[t], {t}];
  DeleteDuplicates[{"exp", NormTerm[#]} & /@ DeleteCases[ts, _?NumericQ]]];
(* Log[2 x] = Log[2] + Log[x], Log[x^3] = 3 Log[x] and Log[x^2+2x+1] =
   2 Log[x+1] are all the generator of their squarefree factors, so a Log
   contributes one key per non-constant irreducible factor of its argument.
   This is exactly what separates Log[x] from Log[1+x] -- the counterexample
   pair above -- so it is load-bearing, not cosmetic. *)
LogKey[a_] := Module[{t = Quiet[Together[a]], ps, fs},
  ps = DeleteCases[Quiet[{Numerator[t], Denominator[t]}], _?NumericQ];
  If[ps === {}, Return[{}]];     (* Log of a constant carries no generator *)
  fs = Flatten[Function[p, If[Quiet[PolynomialQ[p, FreeSyms[p]]] && FreeSyms[p] =!= {},
         Quiet[FactorsOf[p][[All, 1]]], {p}]] /@ ps, 1];
  Sort[DeleteDuplicates[NormPoly /@ DeleteCases[fs, _?NumericQ]]]];

FamHeads = Join[GenHeads, TrigHeads, HypHeads];

(* the family keys of ONE kernel subexpression *)
KernelKeys[e_, x_] := Module[{h = Head[e]},
  Which[
    MatchQ[e, Power[E, _]], ExpKeys[e[[2]]],
    h === PMExp && Length[e] == 1, ExpKeys[e[[1]]],
    Length[e] != 1, {},
    MemberQ[Join[TrigHeads, {Tan, Cot}], h], {{"trig", NormTerm[Together[e[[1]]]]}},
    MemberQ[Join[HypHeads, {Tanh, Coth}], h], {{"trigh", NormTerm[Together[e[[1]]]]}},
    h === Log, {"log", #} & /@ LogKey[e[[1]]],
    h === ArcTan || h === ArcCot, {{"atan", NormTerm[Together[e[[1]]]]}},
    h === ArcTanh || h === ArcCoth, {{"atanh", NormTerm[Together[e[[1]]]]}},
    KeyExistsQ[InvRadical, h], {{SymbolName[h], Together[e[[1]]]}},
    h === ProductLog, {{"W", Together[e[[1]]]}},
    True, {}]];

(* the family keys of a surface expression.  Sqrt[Log[x]] carries Log[x]'s
   key, because Lemma 3.2's flattening gives it the generator u^2 = Log[x];
   a radical of a rational function is its own family.
   The PMExp rewrite is applied HERE rather than relied on from the caller:
   without it 2^x scans as having no family at all, i.e. as rational in x,
   which would absorb it into another block as though it were inert. *)
SurfaceFamilies[u0_, x_] := Module[{u = PMExpHold[u0, x], ks},
  ks = Cases[{u}, e : ((h_[_] /; MemberQ[FamHeads, h]) | Power[E, _] | PMExp[_]) :> KernelKeys[e, x], Infinity];
  ks = Join[ks, Cases[{u}, Power[b_, r_Rational] /; ! IntegerQ[r] :>
         If[Quiet[RationalFunctionQ[b, {x}]], {{"rad", NormPoly[b]}}, SurfaceFamilies[b, x]], Infinity]];
  (* Flatten, not Union @@ -- an empty Union@@{} / Join@@{} stays unevaluated
     and poisons every comparison downstream (measured on 2^x, whose Log[2]
     exponent factor yields no key at all) *)
  Sort[DeleteDuplicates[Flatten[ks, 1]]]];

(* coupling-capable: the block carries a primitive whose derivative lies in
   the shared base C(x) -- a Log or an InvRational head at an argument
   rational in x alone.  Only such blocks can take part in the cross-block
   coupling of the Lemma, so only they need the guard. *)
BasePrimitiveQ[u_, x_] := Cases[{u},
  (h_[a_] /; (h === Log || KeyExistsQ[InvRational, h]) && Quiet[RationalFunctionQ[a, {x}]]),
  Infinity] =!= {};

(* one merge step: fuse the first pair of blocks that share a family *)
BlockMerge1[st_] := Catch[Module[{nb = st[[1]], nk = st[[2]], L},
  L = Length[nb];
  Do[Do[If[Intersection[nk[[i]], nk[[j]]] =!= {},
      Throw[{Append[Delete[nb, {{i}, {j}}], Join[nb[[i]], nb[[j]]]],
             Append[Delete[nk, {{i}, {j}}], Union[nk[[i]], nk[[j]]]]}]],
    {j, i + 1, L}], {i, L}];
  st]];

(* AdditiveBlocks[f, x] -> {{block, keys}, ...}
   The additive terms of f grouped into the connected components of the
   "shares a generator family" graph.  Terms with no family at all are
   rational in x alone; they are absorbed into the cheapest non-rational
   block rather than given one of their own, because a separate rational
   block pays the whole per-call floor for work try_rational does in
   microseconds, and because keeping them together reproduces today's answer
   byte for byte on the corpus cases where the rational part cancels against
   an elementary piece (#134, #199, #204). *)
AdditiveBlocks[f_, x_] := Module[{fe, terms, ks, st, blocks, keys, ratT, idx},
  fe = Expand[f];
  (* An integrand that is not a sum cannot be split, so return before the
     family scan rather than after it.  This keeps the decomposition off the
     critical path of the ~200 corpus cases that are not sums: with the scan
     running unconditionally the median case paid 3-10% for a partition it
     could never use, which showed up as the median rising 0.061 -> 0.063 s
     and "fastest port" dropping from 88 cases to 74 while the heavy tail
     improved. *)
  If[Head[fe] =!= Plus, Return[{{f, {}}}]];
  fe = Expand[PMExpHold[fe, x]];
  If[Head[fe] =!= Plus, Return[{{f, {}}}]];
  terms = List @@ fe;
  ks = SurfaceFamilies[#, x] & /@ terms;
  ratT = Flatten[Position[ks, {}]];
  idx = Complement[Range[Length[terms]], ratT];
  If[idx === {}, Return[{{Total[terms], {}}}]];
  st = FixedPoint[BlockMerge1, {List /@ idx, ks[[idx]]}];
  {blocks, keys} = st;
  If[ratT =!= {},
    With[{cheap = First[Ordering[(Length[keys[[#]]] * 10^6 + LeafCount[Total[terms[[blocks[[#]]]]]]) & /@ Range[Length[blocks]]]]},
      blocks[[cheap]] = Join[blocks[[cheap]], ratT]]];
  Table[{Total[terms[[Sort[blocks[[k]]]]]], keys[[k]]}, {k, Length[blocks]}]];

(* the heads a block may introduce that make its answer more than elementary *)
SpecialHeadsSP = {ExpIntegralEi, LogIntegral, SinIntegral, CosIntegral,
  SinhIntegral, CoshIntegral, Erf, Erfc, Erfi, Gamma, EllipticF, EllipticE, EllipticPi};

(* the {elementary+special, remainder} halves of a surface answer, split on
   the top-level sum exactly as IntegrateSurfacePartial assembles it *)
SplitParts[a_, x_] := Module[{ts = If[Head[a] === Plus, List @@ a, {a}]},
  {Total[DeleteCases[ts, Inactive[Integrate][_, x]]],
   Total[Cases[ts, Inactive[Integrate][r_, x] :> r]]}];
(* a block CLOSED: no remainder and no special function.  By corollary C2 of
   the coupling Lemma such a block cannot take part in a cross-block
   coupling, so it never needs the guard. *)
ClosedElemQ[p_] := p[[2]] === 0 && FreeQ[p[[1]], Alternatives @@ SpecialHeadsSP];

(* SplitIntegratePartial: the additive fast path of IntegrateSurfacePartial.
   Returns {answer, ok} on success, or None to fall through to the joint path.

   The guard is the coupling Lemma above: let N be the blocks that did not
   close to a plain elementary expression and P the coupling-capable ones.
   With |N and P| <= 1 no cross-block rescue can exist and the split is
   complete; with two or more we decline and let the joint path run, because
   that is exactly the Log[1+x]/x + Log[x]/(1+x) shape whose halves are each
   out of class while the sum is Log[x] Log[1+x].

   The merged `ok` is the conjunction of the blocks' own flags -- not a new
   check.  Each block already verified its own answer against its own
   integrand, and D is linear, so a sum of verified pieces is a verified sum.

   The budget is divided among the blocks so that k blocks cannot cost k
   times $SpecialTimeBudget (an outer TimeConstrained cannot interrupt an
   inner one, so the outer wrapper is not a backstop here). *)
SplitIntegratePartial[blocks_, integrandN_, x_, special_, verbose_, samples0_] :=
 Module[{samples = samples0, x0, k = Length[blocks], rs, ps, Nn, Pp, bad, anti, rem},
  If[samples === None,
    x0 = SamplePoint[integrandN, x];
    samples = {If[x0 =!= None, x0, 7/5], 9/4, 13/5}];
  rs = Block[{$SpecialTimeBudget = $SpecialTimeBudget/k},
    Function[b, IntegrateSurfacePartial[b[[1]], x, "Special" -> special,
      "Verbose" -> verbose, "Samples" -> samples, "Additive" -> False]] /@ blocks];
  (* any block that did not return a {answer, flag} pair -- a failure or a
     status -- sends the whole call to the joint path, so a block-level
     "not elementary"/"not in class" is never mistaken for the sum's *)
  If[! AllTrue[rs, ListQ[#] && Length[#] == 2 && BooleanQ[#[[2]]] &], Return[None]];
  ps = SplitParts[#[[1]], x] & /@ rs;
  Nn = Flatten[Position[ClosedElemQ /@ ps, False]];
  Pp = Flatten[Position[BasePrimitiveQ[#[[1]], x] & /@ blocks, True]];
  bad = Intersection[Nn, Pp];
  If[Length[bad] >= 2,
    If[verbose, Print["  additive split declined: ", Length[bad],
      " coupling-capable blocks did not close (cross-block cancellation possible)"]];
    Return[None]];
  anti = Total[ps[[All, 1]]];
  rem = Can[Total[ps[[All, 2]]]];
  If[verbose, Print["  additive split: ", k, " blocks, ",
    Length[Nn], " not closed, ", Length[Pp], " coupling-capable"]];
  (* a COMPLETE answer carries no Inactive[Integrate] at all.  The joint path
     reaches its Inactive[Integrate] line only for a genuine PartialResult and
     returns complete answers earlier through PresentSpecial, so emitting
     Inactive[Integrate][0, x] here would make a complete split answer read as
     partial -- which the corpus judge scores as a different outcome. *)
  {anti + If[rem === 0, 0, Inactive[Integrate][rem, x]], AllTrue[rs, #[[2]] &]}];

(* SplitIntegrateSpecial: the same fast path for the complete-answer entry.
   Returns {answer, ok} on success, or None to fall through to the joint path.

   Complete mode is stricter than partial mode in one way that matters: a
   block may come back as a STATUS ("not in class", "not elementary",
   "failed", or a strict downgrade of an uncertified special answer) rather
   than an answer, and a block's status is NOT the sum's -- the sum may well
   be elementary when a block is not (that is the counterexample).  So any
   non-pair from any block sends the whole call to the joint path, which is
   also what keeps the `nic` cases bit-identical. *)
SplitIntegrateSpecial[blocks_, integrandN_, x_, strict_, verbose_, samples0_] :=
 Module[{samples = samples0, x0, k = Length[blocks], rs, Nn, Pp, bad},
  If[samples === None,
    x0 = SamplePoint[integrandN, x];
    samples = {If[x0 =!= None, x0, 7/5], 9/4, 13/5}];
  rs = Block[{$SpecialTimeBudget = $SpecialTimeBudget/k},
    Function[b, IntegrateSurfaceSpecial[b[[1]], x, "Verbose" -> verbose,
      "Samples" -> samples, "Strict" -> strict, "Additive" -> False]] /@ blocks];
  If[! AllTrue[rs, ListQ[#] && Length[#] == 2 && BooleanQ[#[[2]]] &], Return[None]];
  Nn = Flatten[Position[ClosedElemQ[{#[[1]], 0}] & /@ rs, False]];
  Pp = Flatten[Position[BasePrimitiveQ[#[[1]], x] & /@ blocks, True]];
  bad = Intersection[Nn, Pp];
  If[Length[bad] >= 2,
    If[verbose, Print["  additive split declined: ", Length[bad],
      " coupling-capable blocks did not close elementarily"]];
    Return[None]];
  If[verbose, Print["  additive split: ", k, " blocks, ", Length[Nn], " with special terms"]];
  {Total[rs[[All, 1]]], AllTrue[rs, #[[2]] &]}];

(* ------------------------------------------------------------------ kernels *)
(* Kernel[kind, col, data]: a kernel column, the element col (Trager
   coordinates) with a known antiderivative.  kind:
     "gamma": s, v, omega, cw, shift, J -- omega^k = cw (-v)^(k(s-1)) xi^(kJ), the
              antiderivative cw^(1/k) e^(-shift) Gamma(s, -v) (s = 0: e^(-shift) Ei(v));
     "ell":   g, q, i, first -- g^i D g / y, antiderivative int g^i dg / y;
     "third": g, q, c, p      -- c(g) D g / (p(g) y), by partial fractions.        *)
MakeKernel[kind_, col_, data_: <||>] := Join[<|"kind" -> kind, "col" -> col|>, data];

KernelRepr[K_] := Switch[K["kind"],
  "gamma", "Kernel(Gamma(s=" <> ToStr[K["s"]] <> ", -v), v = " <> ToStr[K["v"]] <> ")",
  "ell", "Kernel(int " <> ToStr[K["g"]] <> "^" <> ToStr[K["i"]] <> " d" <> ToStr[K["g"]] <> "/sqrt(" <> ToStr[K["q"]] <> "))",
  _, "Kernel(int (" <> ToStr[K["c"]] <> ") d" <> ToStr[K["g"]] <> "/((" <> ToStr[K["p"]] <> ") sqrt(" <> ToStr[K["q"]] <> ")))"];

KernelAntiderivative[K_] := Module[{s, k, C, h, A},
  Switch[K["kind"],
    "gamma",
      s = K["s"];
      k = If[s =!= 0, Denominator[s], 1];
      If[s === 0, Return[Exp[-K["shift"]] K["cw"] ExpIntegralEi[K["v"]], Module]];
      C = K["cw"]^(1/k) Exp[-K["shift"]];
      h = If[k == 2, PerfectRoot[-K["v"], k], None];
      If[h =!= None, Return[C Sqrt[Pi] Erfc[h], Module]];        (* Gamma(1/2, h^2) = Sqrt[Pi] Erfc[h]: a smooth branch *)
      C Gamma[s, -K["v"]],
    "ell", FirstKindBack[EllInt[K["g"], K["q"], K["i"]], K["g"], K["q"], Lookup[K, "first", None]],
    _, FirstKindBack[ThirdInt[K["g"], K["q"], K["c"], K["p"]], K["g"], K["q"], Lookup[K, "first", None]]]];

(* PerfectRoot[z, k]: h with h^k = z when z is a constant times a k-th power of a
   rational function (the constant's principal root), else None *)
PerfectRoot[z0_, k_] := Catch[Module[{z = Factor[Cancel[z0]], syms, n, d, h = 1},
  syms = FreeSyms[z];
  If[syms === {}, Throw[None, "pr"]];
  {n, d} = {Numerator[z], Denominator[z]};
  Do[With[{part = pr[[1]], sgn = pr[[2]]},
      If[FreeSyms[part] === {}, h *= (part^(1/k))^sgn; Continue[]];
      Do[If[FreeSyms[fac[[1]]] === {}, h *= (fac[[1]]^fac[[2]])^(sgn/k); Continue[]];
        If[Mod[fac[[2]], k] != 0, Throw[None, "pr"]];
        h *= fac[[1]]^(sgn Quotient[fac[[2]], k]), {fac, FactorList[part]}]],
    {pr, {{n, 1}, {d, -1}}}];
  h], "pr"];

(* the special answer and the partial answer *)
SpecialResult[elementary_, terms_] := <|"type" -> "special", "elementary" -> elementary, "terms" -> terms,
  "certified" -> False, "necessity" -> None, "certificate" -> None, "independence" -> None, "base" -> None|>;
SpecialResultQ[r_] := AssociationQ[r] && Lookup[r, "type", None] === "special";
PartialResult[elementary_, terms_, remainder_] := <|"type" -> "partial", "elementary" -> elementary, "terms" -> terms,
  "remainder" -> remainder, "remainderCertified" -> None|>;
PartialResultQ[r_] := AssociationQ[r] && Lookup[r, "type", None] === "partial";
ResultExpr[r_] := r["elementary"] + Total[(#[[1]] KernelAntiderivative[#[[2]]]) & /@ r["terms"]];
(* an elementary answer: an expression (not a status list, not a result Association) *)
ExprQ[r_] := ! ListQ[r] && ! AssociationQ[r] && r =!= $Failed && r =!= None;

(* ---------------------------------------------------------- exponential data (S1) *)
(* ScalarOrNone[T, u]: u as a scalar expression when its radical coordinates vanish *)
ScalarOrNone[T_, u0_] := With[{u = TPad[T, u0]}, If[AllTrue[Rest[u], Can[#] === 0 &], Can[u[[1]]], None]];

(* VerifySource[T, xi, lam]: D xi / xi == D lambda, exactly *)
VerifySource[T_, xi_, lam_] := QuietCheck[
  With[{lhs = TPad[T, TDiv[T, TowerD[T, xi], xi]], rhs = TPad[T, TowerD[T, lam]]},
    AllTrue[lhs - rhs, Can[#] === 0 &]], False];

(* ExpSources[T, back, Y, verbose]: Algorithm S1, the exponential sources
   {name, xi, lambda} of the tower (xi an exponential element with D xi / xi =
   D lambda) from the surface form (BuildTower's back rules) when available and
   otherwise by inference from the derivation; each verified exactly.  xi =
   e^lambda, exact for the surface generators Exp, Log, ProductLog, Tan, Tanh. *)
ExpSources[T_, back_, Ysym_, verbose_: False] := Module[{gens = T["gens"], q = T["q"], out = {}, flat = <||>, conv,
    old, new, base, mm, olds, s0, symEl, res, xi, lam, nm},
  If[back =!= None,
    (* flattened generators: old symbol -> expression in the current generators *)
    Do[{old, new} = List @@ rule;
      If[MatchQ[new, Power[_, _Rational]] && MemberQ[gens, old],
        {base, mm} = {new[[1]], Denominator[new[[2]]]};
        olds = Select[FreeSyms[base], ! MemberQ[gens, #] && # =!= Ysym &];
        If[Length[olds] == 1 && PolynomialQ[base, olds[[1]]] && Exponent[base, olds[[1]]] == 1,
          s0 = olds[[1]];
          flat[s0] = (old^mm - (base /. s0 -> 0))/Coefficient[base, s0, 1]]],
      {rule, back}];
    conv[e_] := Module[{r = e},
      Do[r = r /. Normal[flat], {Length[flat] + 1}];
      If[q =!= None, r = r /. Power[b_, 1/2] /; Expand[b - q] === 0 :> Ysym];
      If[Ysym =!= None, FromY[T, r, Ysym], TScalar[T, r]]];
    Do[{old, new} = List @@ rule;
      If[old === Ysym, Continue[]];
      symEl = If[MemberQ[gens, old] || KeyExistsQ[flat, old], conv[old], None];
      If[symEl === None, Continue[]];
      res = QuietCheck[Which[
          MatchQ[new, Power[E, _]], {symEl, conv[new[[2]]], "exp(" <> ToStr[new[[2]]] <> ")"},
          Head[new] === Log && Length[new] == 1, {conv[new[[1]]], symEl, ToStr[new[[1]]]},
          Head[new] === ProductLog && Length[new] == 1, {TDiv[T, conv[new[[1]]], symEl], symEl, "exp(W(" <> ToStr[new[[1]]] <> "))"},
          Head[new] === Tan, {TDiv[T, Padd[TScalar[T, 1], Pscale[I, symEl]], Padd[TScalar[T, 1], Pscale[-I, symEl]]],
                              Pscale[2 I, conv[new[[1]]]], "exp(2i " <> ToStr[new[[1]]] <> ")"},
          Head[new] === Tanh, {TDiv[T, Padd[TScalar[T, 1], symEl], Padd[TScalar[T, 1], Pscale[-1, symEl]]],
                               Pscale[2, conv[new[[1]]]], "exp(2 " <> ToStr[new[[1]]] <> ")"},
          True, None], $Failed];
      If[res === None || res === $Failed, Continue[]];
      {xi, lam, nm} = res;
      xi = CanT[TPad[T, xi]]; lam = CanT[TPad[T, lam]];
      If[Ysym =!= None, nm = StringReplace[nm, ToStr[Ysym] -> "y"]];
      If[VerifySource[T, xi, lam],
        AppendTo[out, {nm, xi, lam}];
        If[verbose, Print["  exponential source ", nm, ": xi = ", ToY[T, xi, Global`y], ", lambda = ", ToY[T, lam, Global`y]]]],
      {rule, back}]];
  If[back === None, out = Join[out, InferSources[T, verbose]]];
  out];

(* InferSources[T, verbose]: the sources of a bare Tower -- a hyperexponential
   generator theta with D theta / theta = D eta (eta by Part II below theta), and
   a primitive generator t with D t = D g / g (g from the logarithmic
   antiderivative below t) *)
InferSources[T_, verbose_: False] := Module[{gens = T["gens"], q = T["q"], out = {}, g, d, below, Tb, w, r, Yd, lam, xi, pieces, xie, ok, cf, rest},
  Do[g = gens[[k]]; d = T["derivs"][[k]];
    If[AnyTrue[Rest[d], Can[#] =!= 0 &], Continue[]];
    below = Take[gens, k - 1];
    If[q =!= None && ! SubsetQ[below, Select[gens, ! FreeQ[q, #] &]], Continue[]];
    Tb = Tower[below, Take[T["derivs"], k - 1], If[q =!= None && SubsetQ[below, Select[gens, ! FreeQ[q, #] &]], q, None], T["m"]];
    w = Can[d[[1]]/g];
    QuietCheck[
      If[FreeQ[w, g],
        r = Block[{$analyses = <||>}, ParallelIntegrateMixed[If[Tb["q"] =!= None, TScalar[T, w], {w, 0}], Tb]];
        If[ExprQ[r] && FreeQ[r, _Log | _ArcTan | _ArcTanh | _RootSum | _ExpIntegralEi],
          Yd = Unique["y"];
          lam = FromY[T, If[q =!= None, r /. Sqrt[q] -> Yd, r], Yd];
          xi = TScalar[T, g];
          If[VerifySource[T, xi, lam], AppendTo[out, {ToStr[g], xi, lam}]]],
        If[FreeQ[d[[1]], g],
          r = Block[{$analyses = <||>}, ParallelIntegrateMixed[TScalar[T, d[[1]]], Tb]];
          If[ExprQ[r],
            pieces = If[Head[Expand[r]] === Plus, List @@ Expand[r], {Expand[r]}];
            xie = 1; ok = True;
            Do[{cf, rest} = If[Head[pc] === Times && NumericQ[First[pc]], {First[pc], Rest[pc]}, {1, pc}];
              If[Head[rest] === Log && IntegerQ[cf], xie *= rest[[1]]^cf, ok = False], {pc, pieces}];
            If[ok,
              Yd = Unique["y"];
              xi = FromY[T, If[q =!= None, xie /. Sqrt[q] -> Yd, xie], Yd];
              lam = TScalar[T, g];
              If[VerifySource[T, xi, lam], AppendTo[out, {"exp(" <> ToStr[g] <> ")", xi, lam}]]]]]], Null],
    {k, 2, Length[gens]}];
  If[verbose, Do[Print["  inferred exponential source ", s[[1]]], {s, out}]];
  out];

(* TauEmax[T, tau]: an exponent bound for the exponential monomials whose
   reductions can appear in the decomposition of the residue tau: its degrees in
   the generators, at least 2 *)
TauEmax[T_, tau_] := Module[{e = 2, n, d},
  Do[{n, d} = {Numerator[#], Denominator[#]} &[Can[c]];
    Do[Do[If[! FreeQ[part, g] && PolynomialQ[part, g], e = Max[e, Exponent[part, g]]], {part, {n, d}}], {g, T["gens"]}],
    {c, tau}];
  e];

(* Tpow[T, u, j]: u^j for an integer j *)
Tpow[T_, u_, j_] := Module[{r = TScalar[T, 1], base = If[j > 0, u, TDiv[T, TScalar[T, 1], u]]},
  Do[r = TMul[T, r, base], {Abs[j]}]; TPad[T, r]];

(* Monomials[T, sources, emax, depth]: the exponential monomials prod xi_i^j_i with
   |j_i| <= emax, at most depth nonzero exponents: {exponent vector, xi, lambda} *)
Monomials[T_, sources_, emax_: 2, depth_: 2] := Module[{n = Length[sources], out = {}, xi, lam, ev},
  Do[Do[xi = TScalar[T, 1]; lam = TZero[T];
      Do[xi = TMul[T, xi, Tpow[T, sources[[idx[[t]], 2]], js[[t]]]];
        lam = Padd[lam, Pscale[js[[t]], sources[[idx[[t]], 3]]]], {t, Length[idx]}];
      ev = ConstantArray[0, n]; ev[[idx]] = js;
      AppendTo[out, {ev, CanT[TPad[T, xi]], CanT[TPad[T, lam]]}],
      {js, Tuples[ConstantArray[DeleteCases[Range[-emax, emax], 0], Length[idx]]]}],
    {idx, Join @@ Table[Subsets[Range[n], {r}], {r, 1, depth}]}];
  out];

(* ------------------------------------------------------ residues at primes (S2) *)
(* LinearPlaces[T, p]: the places over the prime p on which the residue is
   evaluated by a substitution: p linear in its main generator g (root rho), or p
   in F[g] with constant coefficients split over the algebraic numbers.
   {{g, rho, ell}, ...} with ell = g - rho, or None *)
LinearPlaces[T_, p_] := Module[{gens = T["gens"], main, g, cl, rts},
  main = Select[gens, ! FreeQ[p, #] &];
  If[main === {}, Return[None]];
  g = Last[main];
  If[! PolynomialQ[p, g], Return[None]];
  If[Exponent[p, g] == 1,
    cl = CoefficientList[p, g];
    Return[{{g, Can[-cl[[1]]/cl[[2]]], g - Can[-cl[[1]]/cl[[2]]]}}]];
  If[FreeQ[CoefficientList[p, g], Alternatives @@ gens],
    rts = RootsMult[p, g];
    If[Total[rts[[All, 2]]] == Exponent[p, g], Return[{g, #, g - #} & /@ rts[[All, 1]]]]];
  None];

(* ResidueAt[T, f, g, rho, ell]: (f ell / D ell) at g = rho, as a tuple in the
   remaining generators (and y) *)
ResidueAt[T_, f_, g_, rho_, ell_] := Module[{lp = TScalar[T, ell], h},
  h = TPad[T, TDiv[T, TMul[T, f, lp], TowerD[T, lp]]];
  Can[Quiet[Can[#] /. g -> rho]] & /@ h];

(* IsConst[T, u, g, rho]: u constant at every place over the prime *)
IsConst[T_, u_, g_: None, rho_: None] := Module[{gens = T["gens"], qv},
  If[! FreeQ[u[[1]], Alternatives @@ gens], Return[False]];
  If[AllTrue[Rest[u], Can[#] === 0 &], Return[True]];
  If[T["q"] === None || T["m"] != 2 || g === None, Return[False]];
  qv = Can[T["q"] /. g -> rho];
  FreeQ[u[[2]], Alternatives @@ gens] && FreeQ[qv, Alternatives @@ gens]];

(* the Ei candidates of a residue: the monomials whose lambda is integral and
   constant at the place, as {ev, xi, lamS, val, K, tK} *)
EiCandidates[T_, monosHere_, g_, rho_, ell_] := Module[{gens = T["gens"], cands = {}, lamS, val, v, vt, K, tK},
  Do[{lamS} = {ScalarOrNone[T, mo[[3]]]};
    If[lamS === None, Continue[]];
    val = Can[Quiet[Can[lamS] /. g -> rho]];
    If[! FiniteValQ[val] || ! FreeQ[val, Alternatives @@ gens], Continue[]];
    v = Can[lamS - val];
    If[v === 0, Continue[]];
    vt = TScalar[T, v];
    K = CanT[TPad[T, TMul[T, TDiv[T, mo[[2]], vt], TowerD[T, vt]]]];       (* xi D v / v *)
    tK = ResidueAt[T, K, g, rho, ell];
    If[AllTrue[tK, Can[#] === 0 &], Continue[]];
    AppendTo[cands, {mo[[1]], mo[[2]], lamS, val, K, tK}],
    {mo, monosHere}];
  cands];

(* DecompResidue[T, tt, cands]: r with tt = r0 + Sum r_j tK_j coordinatewise in
   the residue field, the free unknowns 0; None if tt is outside the span *)
DecompResidue[T_, tt_, cands_] := Module[{gens = T["gens"], r0 = Unique["r0"], rs, eqs = {}},
  rs = Table[Unique["r"], {Length[cands]}];
  Do[eqs = Join[eqs, NumCoeffs[tt[[i]] - If[i == 1, r0, 0] - Sum[rs[[j]] cands[[j, 6, i]], {j, Length[cands]}], gens]],
    {i, Length[tt]}];
  If[eqs === {}, None, LinSolveZero[eqs, Prepend[rs, r0]]]];

(* EiFromResidues[f, T, sources, verbose]: Algorithm S2, the s = 0 kernels
   determined by residues (Theorem 4.2): {f remaining, {{coefficient, K}, ...},
   certificate or None} *)
EiFromResidues[f0_, T_, sources_, verbose_: False] := Catch[Module[{f = f0, gens = T["gens"], q = T["q"], nc, den = 1,
    monos, found = {}, p, branch, eta, delta, special, vpf, places, g, rho, ell, tr, tau, em, monosHere, cands, sub, ev, xi, lamS, val, K, v, kern, cert},
  nc = NCo[T];
  Do[den = PolynomialLCM[den, Denominator[Can[f[[i]]]]], {i, nc}];
  If[FreeQ[den, Alternatives @@ gens], Throw[{f, {}, None}, "efr"]];
  monos = If[sources =!= {}, Monomials[T, sources], {}];
  Do[p = fac[[1]];
    If[FreeQ[p, Alternatives @@ gens], Continue[]];
    {branch, eta, delta, special} = ClassifyPrime[T, p];
    If[special || branch, Continue[]];
    vpf = VP[f, p, T, branch];
    If[vpf > -delta, Continue[]];
    places = LinearPlaces[T, p];
    If[places === None, Continue[]];
    If[vpf < -delta && (q =!= None || Length[places] != 1 || Expand[places[[1, 3]] - p] =!= 0), Continue[]];   (* deep poles: transcendental towers, linear primes *)
    Do[{g, rho, ell} = pl;
      If[vpf < -delta,
        tr = CanonicalResidueField[T, f, p, g, delta, False];          (* Algorithm 2 of Part II (local Hermite reduction) *)
        If[tr === None, Continue[]];
        tau = TScalar[T, Can[tr]],
        tau = ResidueAt[T, f, g, rho, ell]];
      If[IsConst[T, tau, g, rho], Continue[]];
      (* candidate monomials: lambda integral and constant at the place; the exponent bound read off tau *)
      em = TauEmax[T, tau];
      monosHere = If[em <= 2, monos, If[sources =!= {}, Monomials[T, sources, em], {}]];
      cands = EiCandidates[T, monosHere, g, rho, ell];
      sub = DecompResidue[T, tau, cands];
      If[sub === None,
        cert = {"not in class", p, tau,
          "canonical residue outside Fbar + span of the reductions of the exponential monomials constant at the prime (Corollary 4.3; relative to the sources " <>
          ToStr[sources[[All, 1]]] <> ")"};
        If[verbose, Print["  (", p, ") at ", g, " = ", rho, ": residue ", tau, " is not in the exponential span"]];
        Throw[{f0, {}, cert}, "efr"]];
      Do[If[ExactZeroQ[sub[[j + 1]]], Continue[]];
        {ev, xi, lamS, val, K} = cands[[j, 1 ;; 5]];
        v = Can[lamS - val];
        (* xi = e^lambda = e^(v + val): the antiderivative e^val Ei(v), stored as shift = -val *)
        kern = MakeKernel["gamma", K, <|"s" -> 0, "v" -> v, "omega" -> CanT[TPad[T, TDiv[T, xi, TScalar[T, v]]]], "cw" -> 1,
          "shift" -> -val, "J" -> ev, "place" -> {g, rho}|>];
        AppendTo[found, {sub[[j + 1]], kern}];
        f = Padd[f, Pscale[-sub[[j + 1]], K]];
        If[verbose, Print["  (", p, ") at ", g, " = ", rho, ": residue ", tau, " -> Ei kernel with v = ", v, ", coefficient ", sub[[j + 1]]]],
        {j, Length[cands]}],
      {pl, places}],
    {fac, FactorsOf[den]}];
  {CanRT[f], found, None}], "efr"];

(* ---------------------------------------------------------- Gamma kernels (S3) *)
(* ExactRoot[W, k, gens]: {R, c} with W = c R^k, c constant, if W is a constant
   times a k-th power of a rational function of gens; else None *)
ExactRoot[W0_, k_, gens_] := Catch[Module[{W = Can[W0], R = 1, c = 1, fl},
  (* FactorList over Q leaves a GAUSSIAN k-th power irreducible -- for the
     Fresnel pair Sin[x^2] / Cos[x^2] the square W = I xi^2/x^2 arrives as
     (I - 2 t - I t^2)/(x^2 - 2 I t x^2 - t^2 x^2), whose numerator is
     -I (t - I)^2 but factors over Q as itself, multiplicity 1.  1 is not
     divisible by k = 2, so the root was declined and the Erf kernel never
     offered.  Factor over Q(i) when the coefficients are not rational.
     Extension -> Automatic is NOT enough here (measured: it hands the
     polynomial back unfactored); the generator has to be named.
     The constant this frees needs no root of its own -- it accumulates in c and
     is returned as the kernel's cw, the identity being omega^k = W/cw -- which
     is why Sqrt[I] never has to enter the constant field. *)
  fl[u_] := If[FreeQ[u, _Complex], FactorList[u], Quiet[FactorList[u, Extension -> I]]];
  Do[With[{part = pr[[1]], sgn = pr[[2]]},
      If[FreeQ[part, Alternatives @@ gens], c *= part^sgn; Continue[]];
      Do[If[FreeQ[fac[[1]], Alternatives @@ gens], c *= (fac[[1]]^fac[[2]])^sgn; Continue[]];
        If[Mod[fac[[2]], k] != 0, Throw[None, "er"]];
        R *= fac[[1]]^(sgn Quotient[fac[[2]], k]), {fac, fl[part]}]],
    {pr, {{Numerator[W], 1}, {Denominator[W], -1}}}];
  {R, c}], "er"];

(* HypSourceQ[T, xi]: the source xi is a hyperexponential generator theta != x
   (D theta = lambda theta with lambda free of theta, no radical coordinate) *)
HypGen[T_, xi_] := Module[{gens = T["gens"], th = ScalarOrNone[T, xi], k},
  If[th === None || ! MemberQ[gens, th] || th === gens[[1]], Return[None]];
  k = Position[gens, th][[1, 1]];
  If[AllTrue[Rest[T["derivs"][[k]]], Can[#] === 0 &] && FreeQ[Can[T["derivs"][[k, 1]]/th], th], th, None]];

(* Support[T, f, sources]: the exponents of each hyperexponential-type source in
   f (a weight filter for the Kummer search), with the fractional weight carried
   by y when the radicand is monomial in the source *)
Support[T_, f_, sources_] := Module[{gens = T["gens"], q = T["q"], sup = <||>, th, js, wy, P, comp, n, d, a, jj},
  Do[th = HypGen[T, sources[[si, 2]]];
    js = {};
    If[th =!= None,
      wy = 0;
      If[q =!= None && ! FreeQ[q, th],
        P = CoefficientRules[q, {th}];
        If[Length[P] == 1, wy = If[T["m"] == 2, Exponent[q, th]/2, None]]];
      Do[comp = Can[f[[i]]];
        If[comp === 0, Continue[]];
        {n, d} = {Numerator[comp], Denominator[comp]};
        If[Length[CoefficientRules[d, {th}]] != 1, js = Union[js, {1, -1}]; Continue[]];
        a = Exponent[d, th];
        Do[jj = mon[[1]] - a + If[wy =!= None, (i - 1) wy, 0];
          If[jj != 0, js = Union[js, {jj}]], {mon, Keys[CoefficientRules[n, {th}]]}],
        {i, NCo[T]}],
      js = {-1, 1}];
    sup[si] = Sort[js],
    {si, Length[sources]}];
  sup];

(* TpowFrac[T, xi, J]: xi^J for an integer J; for a half-integer J only when xi
   is a generator theta and y^2 = q is monomial in theta: theta^J sqrt(q0) as
   y theta^(J - a/2), q = c theta^a *)
TpowFrac[T_, xi_, J_] := Module[{th, q = T["q"], a, e},
  If[IntegerQ[J], Return[Tpow[T, xi, J]]];
  th = ScalarOrNone[T, xi];
  If[th === None || q === None || T["m"] != 2 || Denominator[J] != 2 || ! MemberQ[T["gens"], th], Return[None]];
  If[Length[CoefficientRules[q, {th}]] != 1 || EvenQ[Exponent[q, th]], Return[None]];
  a = Exponent[q, th];
  e = J - a/2;
  If[! IntegerQ[e], Return[None]];
  TPad[T, TMul[T, {0, 1}, Tpow[T, TScalar[T, th], e]]]];

(* TowerSpecialPlaces[T]: {g, rho} for the special primes of the tower that are
   linear in their main generator: the hyperexponential and Lambert generators,
   the arguments of logarithms, found as the special factors of the derivatives'
   numerators and denominators and of the generators *)
TowerSpecialPlaces[T_] := Module[{gens = T["gens"], cand, out = {}, pl},
  cand = Rest[gens];
  Do[Do[If[Can[comp] === 0, Continue[]];
      Do[If[FreeQ[part, Alternatives @@ gens], Continue[]];
        cand = Join[cand, FactorsOf[part][[All, 1]]], {part, {Numerator[Can[comp]], Denominator[Can[comp]]}}],
      {comp, d}], {d, T["derivs"]}];
  Do[If[! TrueQ[QuietCheck[ClassifyPrime[T, pp][[4]], False]], Continue[]];
    pl = LinearPlaces[T, pp];
    If[pl =!= None && pl =!= {}, out = Join[out, pl[[All, 1 ;; 2]]]],
    {pp, DeleteDuplicates[cand]}];
  out];

(* JointMonomials[T, f, sources]: the exponential monomials prod theta_k^j_k with
   at least two hyperexponential sources (or one with a logarithmic source) that
   occur in f, as combined sources {name, xi, lambda} *)
JointMonomials[T_, f_, sources_] := Module[{gens = T["gens"], hyp = {}, logsrc, ths, out = {}, seen = {}, comp, n, d, Pd, dmon, J,
    lsym, dl, key, xi, lam, th},
  Do[th = HypGen[T, src[[2]]];
    If[th =!= None, AppendTo[hyp, {th, src[[2]], src[[3]]}]], {src, sources}];
  (* the arguments of logarithmic generators are exponential sources of weight 0 for the
     hyperexponential grading: e^(-log(x)^2) x = e^(-(t - 1/2)^2 + 1/4) *)
  logsrc = Select[sources, Function[src, ! AnyTrue[hyp, #[[2]] === src[[2]] &] && ScalarOrNone[T, src[[3]]] =!= None &&
      MemberQ[gens, ScalarOrNone[T, src[[3]]]]]];
  ths = hyp[[All, 1]];
  If[hyp =!= {} && logsrc =!= {},
    Do[comp = Can[f[[i]]];
      If[comp === 0, Continue[]];
      {n, d} = {Numerator[comp], Denominator[comp]};
      Pd = CoefficientRules[d, ths];
      If[Length[Pd] != 1, Continue[]];
      dmon = Pd[[1, 1]];
      Do[J = mon - dmon;
        If[! AnyTrue[J, # != 0 &], Continue[]];
        Do[lsym = ScalarOrNone[T, ls[[2]]];
          dl = Max[Join[{0}, If[MemberQ[gens, lsym], Exponent[#, lsym] & /@ Select[{n, d}, ! FreeQ[#, lsym] &], {}]]];
          Do[key = Join[J, {ls[[1]], jl}];
            If[MemberQ[seen, key], Continue[]];
            AppendTo[seen, key];
            xi = Tpow[T, ls[[2]], jl]; lam = Pscale[jl, ls[[3]]];
            Do[If[J[[t]] != 0,
                xi = TMul[T, xi, Tpow[T, hyp[[t, 2]], J[[t]]]];
                lam = Padd[lam, Pscale[J[[t]], hyp[[t, 3]]]]], {t, Length[hyp]}];
            AppendTo[out, {"prod " <> ToStr[J] <> " x " <> ls[[1]] <> "^" <> ToString[jl], CanT[TPad[T, xi]], CanT[TPad[T, lam]]}],
            {jl, DeleteCases[Range[-(dl + 2), dl + 2], 0]}],
          {ls, logsrc}],
        {mon, Keys[CoefficientRules[n, ths]]}],
      {i, NCo[T]}]];
  If[Length[hyp] < 2, Return[out]];
  Do[comp = Can[f[[i]]];
    If[comp === 0, Continue[]];
    {n, d} = {Numerator[comp], Denominator[comp]};
    Pd = CoefficientRules[d, ths];
    If[Length[Pd] != 1, Continue[]];
    dmon = Pd[[1, 1]];
    Do[J = mon - dmon;
      If[Count[J, _?(# != 0 &)] < 2 || MemberQ[seen, J], Continue[]];
      AppendTo[seen, J];
      xi = TScalar[T, 1]; lam = TZero[T];
      Do[If[J[[t]] != 0,
          xi = TMul[T, xi, Tpow[T, hyp[[t, 2]], J[[t]]]];
          lam = Padd[lam, Pscale[J[[t]], hyp[[t, 3]]]]], {t, Length[hyp]}];
      AppendTo[out, {"prod " <> ToStr[J], CanT[TPad[T, xi]], CanT[TPad[T, lam]]}],
      {mon, Keys[CoefficientRules[n, ths]]}],
    {i, NCo[T]}];
  out];

(* GammaCandidates[f, T, sources, verbose]: Algorithm S3, s > 0 (and s = 0 at
   silent zeros): the Kummer search of Section 6.4 of the paper, for arguments
   that are rational functions of a single generator *)
GammaCandidates[f_, T_, sources0_, verbose_: False] := Module[{gens = T["gens"], q = T["q"], m = T["m"], out = {}, seen = {},
    sup, specialPlaces, sources = sources0, lamS, Lam, vars, g, N0, Dn, C, shifts, silent, PN, res, val, lim, v, xiJ, key, vt, K,
    num, mults, Nc, ks, s, kJ, W, Wxi, WxiS, cands, rt0, rt1, omega, cw, lhs, rhs, allc, inSetQ},
  inSetQ[c_, set_] := AnyTrue[set, ExactZeroQ[# - c] &];
  sup = Support[T, f, sources];
  specialPlaces = TowerSpecialPlaces[T];
  Do[AppendTo[sources, comb]; sup[Length[sources]] = {1}, {comb, JointMonomials[T, f, sources]}];   (* products of hyperexponentials: e^(x^2) e^x *)
  Do[lamS = ScalarOrNone[T, sources[[si, 3]]];
    If[lamS === None, Continue[]];
    Do[Lam = Can[J lamS];
      vars = Select[gens, ! FreeQ[Lam, #] &];
      If[Length[vars] != 1, Continue[]];
      g = vars[[1]];
      {N0, Dn} = {Numerator[Together[Lam]], Denominator[Together[Lam]]};
      C = Unique["C"];
      shifts = {};
      QuietCheck[
        PN = Expand[N0 + C Dn];
        If[PolynomialQ[PN, g] && Exponent[PN, g] >= 1,
          res = Resultant[PN, D[PN, g], g];
          If[! FreeQ[res, C],
            Do[If[! inSetQ[r, shifts], AppendTo[shifts, r]], {r, RootsOf[res, C]}]]], Null];
      If[q =!= None && ! FreeQ[q, g],
        Do[If[! FreeQ[rq, Alternatives @@ gens], Continue[]];
          val = Quiet[Together[Lam /. g -> rq]];
          If[FiniteValQ[val] && ! inSetQ[-val, shifts], AppendTo[shifts, -val]],
          {rq, RootsOf[q, g]}]];
      (* silent zeros for s = 0: special generators and infinity *)
      silent = {};
      Do[If[! FreeQ[Lam, sp[[1]]],
          val = QuietCheck[Can[Can[Lam] /. sp[[1]] -> sp[[2]]], $Failed];
          If[val =!= $Failed && FreeQ[val, Alternatives @@ gens] && FiniteValQ[val] && ! inSetQ[-val, silent], AppendTo[silent, -val]]],
        {sp, specialPlaces}];
      lim = QuietCheck[Limit[Lam, g -> Infinity], $Failed];
      If[lim =!= $Failed && FiniteValQ[lim] && FreeQ[lim, Alternatives @@ gens] && FreeQ[lim, Limit] && ! inSetQ[-lim, silent], AppendTo[silent, -lim]];
      allc = shifts;
      Do[If[! inSetQ[c, allc], AppendTo[allc, c]], {c, silent}];
      Do[v = Can[Lam + c];
        If[v === 0, Continue[]];
        xiJ = TpowFrac[T, sources[[si, 2]], J];
        (* s = 0 at silent zeros: xi^J D v / v *)
        If[inSetQ[c, silent] && xiJ =!= None && IntegerQ[J],
          key = {"ei", v};
          If[! MemberQ[seen, key],
            AppendTo[seen, key];
            vt = TScalar[T, v];
            K = CanT[TPad[T, TMul[T, TDiv[T, xiJ, vt], TowerD[T, vt]]]];
            AppendTo[out, MakeKernel["gamma", K, <|"s" -> 0, "v" -> v, "omega" -> CanT[TPad[T, TDiv[T, xiJ, vt]]], "cw" -> 1,
              "shift" -> c, "J" -> J|>]];
            If[verbose, Print["  silent Ei kernel: v = ", v]]]];
        If[! inSetQ[c, shifts], Continue[]];
        (* the multiplicities of the zeros of v *)
        num = Numerator[Together[v]];
        mults = {};
        Do[With[{br = q =!= None && ! FreeQ[q, g] && ExactZeroQ[Quiet[Together[q /. g -> rm[[1]]]]]},
            AppendTo[mults, rm[[2]] If[br, m, 1]]],
          {rm, RootsMult[num, g]}];
        Nc = If[mults === {}, 0, GCD @@ mults];
        ks = If[Nc != 0, Select[Divisors[Nc], # >= 2 &], {}];
        If[ks === {}, ks = If[q =!= None, {2}, {}]];                   (* Cor. 3.2: k | N_c *)
        Do[Do[If[GCD[a, k] != 1, Continue[]];
            s = a/k; kJ = k J;
            If[! IntegerQ[kJ], Continue[]];
            W = (-v)^(k (s - 1));
            Wxi = TpowFrac[T, sources[[si, 2]], kJ];
            If[Wxi === None, Continue[]];
            WxiS = ScalarOrNone[T, Wxi];
            If[WxiS === None, Continue[]];
            W = Can[W WxiS];
            cands = {};
            rt0 = ExactRoot[W, k, gens];
            If[rt0 =!= None, AppendTo[cands, {TScalar[T, rt0[[1]]], rt0[[2]]}]];
            If[q =!= None && m == 2 && k == 2,
              rt1 = ExactRoot[Can[W/q], 2, gens];
              If[rt1 =!= None, AppendTo[cands, {TPad[T, {0, rt1[[1]]}], rt1[[2]]}]]];
            Do[{omega, cw} = cd; omega = CanT[omega];
              vt = TScalar[T, v];
              lhs = TPad[T, TDiv[T, TowerD[T, omega], omega]];
              rhs = TPad[T, Padd[TowerD[T, vt], Pscale[s - 1, TDiv[T, TowerD[T, vt], vt]]]];
              If[! AllTrue[lhs - rhs, Can[#] === 0 &], Continue[]];
              key = {"g", s, v, omega};
              If[MemberQ[seen, key], Continue[]];
              AppendTo[seen, key];
              (* omega^k = cw' W with W = (-v)^(k(s-1)) xi^(kJ): cw' = 1/cw *)
              K = CanT[TPad[T, TMul[T, omega, TowerD[T, vt]]]];
              AppendTo[out, MakeKernel["gamma", K, <|"s" -> s, "v" -> v, "omega" -> omega, "cw" -> 1/cw, "shift" -> c, "J" -> J|>]];
              If[verbose, Print["  Gamma kernel: s = ", s, ", v = ", v, ", omega = ", ToY[T, omega, Global`y]]],
              {cd, cands}],
            {a, 1, k - 1}],
          {k, ks}],
        {c, SortBy[allc, {N[Re[#]], N[Im[#]], LeafCount[#]} &]}],
      {J, Lookup[sup, si, {}]}],
    {si, Length[sources]}];
  out];

(* -------------------------------------------------------- elliptic kernels (S4) *)
(* Pencil[T]: {g, r} when y^2 = q(g) with constant coefficients, deg >= 3, and
   D g = r(g) with no radical coordinate: the pencil of Section 7 *)
Pencil[T_] := Module[{q = T["q"], gens = T["gens"], d},
  If[q === None || T["m"] != 2, Return[None]];
  Do[If[SubsetQ[{gens[[k]]}, FreeSyms[q]] && PolynomialQ[q, gens[[k]]] && Exponent[q, gens[[k]]] >= 3 &&
        FreeSyms[CoefficientList[q, gens[[k]]]] === {},
      d = T["derivs"][[k]];
      If[AllTrue[Rest[d], Can[#] === 0 &], Return[{gens[[k]], Can[d[[1]]]}, Module]]],
    {k, Length[gens]}];
  None];

EllipticColumns[T_, verbose_: False] := Module[{pen = Pencil[T], g, r, q = T["q"], d, first, out = {}},
  If[pen === None, Return[{}]];
  {g, r} = pen;
  d = Exponent[q, g];
  first = FirstKindGenerator[T, g, r];
  Do[AppendTo[out, MakeKernel["ell", TPad[T, {0, Can[g^i r/q]}], <|"g" -> g, "q" -> q, "i" -> i, "first" -> first|>]],   (* g^i D g / y = g^i r y / q *)
    {i, 0, d - 2}];
  If[verbose, Print["  elliptic pencil ", g, ": y^2 = ", q, ", ", d - 1, " columns g^i Dg/y"]];
  out];

(* FirstKindGenerator[T, g, r]: {z, kappa} when a generator z of the tower has
   D z = kappa D g / y: the first-kind integral is then the element z / kappa *)
FirstKindGenerator[T_, g_, r_] := Module[{col1 = Can[r/T["q"]], dz, kap},
  Do[dz = T["derivs"][[k]];
    If[T["gens"][[k]] === g || Can[dz[[1]]] =!= 0, Continue[]];
    kap = Can[dz[[2]]/col1];
    If[kap =!= 0 && FreeSyms[kap] === {}, Return[{T["gens"][[k]], kap}, Module]],
    {k, Length[T["gens"]]}];
  None];

(* ThirdKind[f, T, verbose]: Algorithm S4, Steps 4--7: odd residues at simple
   poles off the branch locus of a pencil curve, removed by c(g) D g / (p y),
   c = f_1 p q / D g mod p *)
ThirdKind[f0_, T_, verbose_: False] := Module[{f = f0, pen = Pencil[T], g, r, q = T["q"], f1, n, d, found = {}, p, h, hn, hd, eg, inv, c, K},
  If[pen === None, Return[{f, {}}]];
  {g, r} = pen;
  f1 = Can[f[[2]]];
  If[f1 === 0 || ! SubsetQ[{g}, FreeSyms[f1]], Return[{f, {}}]];
  {n, d} = {Numerator[f1], Denominator[f1]};
  Do[p = fac[[1]];
    If[FreeQ[p, g] || fac[[2]] != 1 || PolynomialRemainder[q, p, g] === 0, Continue[]];
    h = Can[f1 p q/r];
    {hn, hd} = {Numerator[h], Denominator[h]};
    If[PolynomialRemainder[hd, p, g] === 0, Continue[]];
    eg = PolynomialExtendedGCD[hd, p, g];
    inv = eg[[2, 1]]/eg[[1]];
    c = Expand[PolynomialRemainder[Expand[hn inv], p, g]];
    If[c === 0, Continue[]];
    K = TPad[T, {0, Can[c r/(p q)]}];
    AppendTo[found, {1, MakeKernel["third", K, <|"g" -> g, "q" -> q, "c" -> c, "p" -> p|>]}];
    f = Padd[f, Pscale[-1, K]];
    If[verbose, Print["  third-kind kernel at (", p, "): c = ", c]],
    {fac, FactorsOf[d]}];
  {CanRT[f], found}];

(* ----------------------------------------------------- Legendre presentation *)
(* FirstKindBack[A, g, q, first]: the first-kind integral int dg / y of the
   presentation written as the tower element z / kappa when the tower contains it *)
FirstKindBack[A_, g_, q_, first_] := Module[{z, kap, J0, Fs, Fa, alpha, beta},
  If[first === None, Return[A]];
  {z, kap} = first;
  J0 = EllInt[g, q, 0];
  Fs = DeleteDuplicates[Cases[{J0}, _EllipticF, Infinity]];
  If[Length[Fs] != 1, Return[A]];
  Fa = Fs[[1]];
  alpha = Coefficient[Expand[J0], Fa];
  beta = Simplify[J0 - alpha Fa];
  If[alpha === 0 || ! FreeQ[beta, Fa], Return[A]];
  Expand[A /. Fa -> (z/kap - beta)/alpha]];

CubicData[g_, q_] := Module[{lc, rts, e1, e2, e3, lam, mm, phi},
  lc = Coefficient[q, g, Exponent[q, g]];
  rts = RootsOf[q, g];
  If[Length[rts] != 3, Return[None]];
  If[AllTrue[rts, Im[N[#]] == 0 &], rts = SortBy[rts, -N[#] &]];
  {e1, e2, e3} = rts;
  lam = Sqrt[e1 - e3];
  mm = Simplify[(e2 - e3)/(e1 - e3)];
  phi = ArcSin[Sqrt[(e1 - e3)/(g - e3)]];
  {lc, e1, e2, e3, lam, mm, phi}];

(* LegendreQuartic[g, q]: {c, k} when q = c (1 - g^2)(1 - k g^2), the Legendre normal form *)
LegendreQuartic[g_, q_] := Module[{c0, Q, kk},
  If[! PolynomialQ[q, g] || Exponent[q, g] != 4, Return[None]];
  c0 = q /. g -> 0;
  If[c0 === 0, Return[None]];
  Q = Expand[q/c0];
  If[Coefficient[Q, g, 1] =!= 0 || Coefficient[Q, g, 3] =!= 0, Return[None]];
  kk = Coefficient[Q, g, 4];
  If[Expand[Q - (1 - g^2) (1 - kk g^2)] =!= 0, Return[None]];
  {c0, kk}];

(* EllInt[g, q, i]: int g^i dg / Sqrt[q]: Legendre form on a cubic for i = 0, 1
   (valid for g beyond the largest real root; the branch fixed at the surface),
   the Legendre quartic for i = 0, 2, otherwise unevaluated *)
EllInt[g_, q_, i_] := Module[{leg, c0, kk, F, E0, cd, lc, e1, e2, e3, lam, mm, phi, y},
  leg = LegendreQuartic[g, q];
  If[leg =!= None && MemberQ[{0, 2}, i],
    {c0, kk} = leg;
    F = EllipticF[ArcSin[g], kk]; E0 = EllipticE[ArcSin[g], kk];
    Return[If[i == 0, F, (F - E0)/kk]/Sqrt[c0]]];
  cd = If[Exponent[q, g] == 3, CubicData[g, q], None];
  If[cd === None || i > 1, Return[Inactive[Integrate][g^i/Sqrt[q], g]]];
  {lc, e1, e2, e3, lam, mm, phi} = cd;
  F = EllipticF[phi, mm]; E0 = EllipticE[phi, mm];
  y = Sqrt[q];
  If[i == 0, -2/lam F/Sqrt[lc],
    (-2/lam (e3 + lam^2) F + 2 lam E0)/Sqrt[lc] + 2 y/(lc (g - e3))]];

(* ThirdInt[g, q, c, p]: int c(g) dg / (p(g) Sqrt[q]) by partial fractions over
   the roots of p: Legendre Pi on a cubic, otherwise unevaluated *)
ThirdInt[g_, q_, c_, p_] := Module[{cd, lc, e1, e2, e3, lam, mm, phi, out = 0, dp, w, n, Pi0},
  cd = If[Exponent[q, g] == 3, CubicData[g, q], None];
  If[cd === None, Return[Inactive[Integrate][c/(p Sqrt[q]), g]]];
  {lc, e1, e2, e3, lam, mm, phi} = cd;
  dp = D[p, g];
  Do[w = Simplify[(c /. g -> rho)/(dp /. g -> rho)];
    n = Simplify[(rho - e3)/lam^2];
    Pi0 = 2/(lam (rho - e3)) (EllipticF[phi, mm] - EllipticPi[n, phi, mm])/Sqrt[lc];
    out += w Pi0,
    {rho, RootsOf[p, g]}];
  out];

(* Present[e]: rewrite into the functions a user expects (additive constants
   dropped, since every term has a constant coefficient): Gamma[1/2, z] through
   Erfc, Erfc[I z] through Erfi, ExpIntegralEi[Log[a]] through LogIntegral, and a
   conjugate pair a Ei(I z) + b Ei(-I z) through CosIntegral and SinIntegral *)
Present[e0_] := Module[{e = e0, imagArg, ats, z, a, b},
  e = e /. Gamma[1/2, z_] :> Sqrt[Pi] Erfc[Sqrt[z]];
  imagArg[w0_] := With[{w = Expand[w0/I]}, If[FreeQ[w, _Complex], w, None]];
  e = e /. Erfc[w_] /; imagArg[w] =!= None :> -I Erfi[imagArg[w]];
  e = e /. ExpIntegralEi[Log[w_]] :> LogIntegral[w];
  e = Expand[e];
  ats = DeleteDuplicates[Cases[{e}, _ExpIntegralEi, Infinity]];
  Do[z = Expand[at[[1]]/I];
    If[! FreeQ[z, _Complex] || MinusSignQ[z], Continue[]];
    a = Coefficient[e, at];
    b = Coefficient[e, ExpIntegralEi[-I z]];
    If[a === 0, Continue[]];
    e = Expand[e - a at - b ExpIntegralEi[-I z] + Simplify[a + b] CosIntegral[z] + Simplify[I (a - b)] SinIntegral[z]],
    {at, ats}];
  e];

(* HalfAngleFold[e, back]: undo the Weierstrass substitution the tower left in
   the surface.  A tangent generator t = Tan[u] makes the answer a rational
   function of Tan[u], so Integrate[Log[x] Sin[x], x] reads

     CosIntegral[x] - Log[x]/(1 + Tan[x/2]^2) + (Log[x] Tan[x/2]^2)/(1 + Tan[x/2]^2)

   which is CosIntegral[x] - Cos[x] Log[x].  Correct, verified, and unreadable.
   (The Python and Maxima references leave it in this form too; the corpus judge
   only asks for a verified answer carrying a special-function head.  This is a
   REPL-quality fix, and it matters because these answers reach plain Integrate.)

   Rewrite EXACTLY the generators named in back -- never an unrelated Tan the
   answer earned honestly, which Sin[2u]/(1 + Cos[2u]) would only make worse --
   clear the fraction, and reduce the even powers of Sin by sin^2 = 1 - cos^2 so
   Cancel can finish.  Deterministic and cheap: no Simplify, which on a
   trigonometric fraction is both slow and a known hang.  Accepted only when it
   strictly shrinks the expression, so the fold can never make output worse; the
   caller verifies the folded form, not the raw one. *)
HalfAngleFold[e_, back_] := Module[{tans, rules, f},
  If[FreeQ[e, Tan], Return[e]];
  tans = DeleteDuplicates[Cases[{back}, _Tan, Infinity]];
  If[tans === {}, Return[e]];
  rules = (# -> Sin[2 First[#]]/(1 + Cos[2 First[#]])) & /@ tans;
  f = QuietCheck[
        Cancel[Together[Together[e /. rules] /. Sin[u_]^n_ /; EvenQ[n] :> (1 - Cos[u]^2)^(n/2)]],
        $Failed];
  (* the test is on the TARGETED generators, not on Tan at large: an unrelated
     Tan the answer earned honestly may survive and the fold still be a win *)
  If[f === $Failed || ! FreeQ[f, Alternatives @@ tans] || LeafCount[f] >= LeafCount[e],
    e, Expand[f]]];

(* ------------------------------------------------------- the linear system (S6) *)
(* SPAnalyse[f, T, verbose]: Part II's analysis (Steps 1--14 of Algorithm 4),
   fresh (no cached analysis of an earlier run is reused, as _analyse), kept in
   $analyses under a key of its own; returns the key, or the status list of an
   early exit *)
SPAnalyse[f_, T_, verbose_: False] := Module[{r, key},
  r = Block[{$analyses = <||>}, With[{k = Analyse[f, T, verbose]}, If[ListQ[k], k, $analyses[k]]]];
  If[ListQ[r], Return[r]];
  key = "sp" <> ToString[++$spCounter];
  $analyses[key] = r;
  key];

(* KernelProxy[key, kernels]: the analysis with the residual replaced by a generic
   combination of the residual and the kernel columns, the bound caches reset *)
KernelProxy[key_, kernels_] := Module[{A = $analyses[key], T, comb, key2},
  T = A["T"]; comb = A["rem"];
  Do[comb = Padd[comb, Pscale[3 + 2 (j - 1), TPad[T, kernels[[j]]["col"]]]], {j, Length[kernels]}];
  key2 = "proxy" <> ToString[++$spCounter];
  $analyses[key2] = Join[KeyDrop[A, {"inf"}], <|"rem" -> CanT[comb], "spec" -> <||>|>];
  key2];

(* ProxyBounds[key, unkLogs, kernels, retry, verbose]: Algorithm 6 with v(f)
   replaced by min(v(rem), v(K)) (Corollary 4.4): the bounds taken on a generic
   combination of the residual and the kernel columns *)
ProxyBounds[key_, unkLogs_, kernels_, retry_, verbose_] := ExtendedBounds[$ExtendedBoundsSP,
  Module[{k2, A2, r},
    k2 = If[kernels === {}, key, KernelProxy[key, kernels]];
    A2 = $analyses[k2];
    r = PlaceBounds[k2, A2["T"], A2["rem"], A2["denv"], unkLogs, retry, verbose];
    If[k2 =!= key, $analyses = KeyDrop[$analyses, {k2}]];
    r]];

(* SolveKernelSystem[T, rem, denv, units, unkLogs, css, monos, gammas, betas, unks, extra]:
   the system of the ansatz with the extra (kernel, remainder) columns: the
   particular solution with the non-pivot unknowns 0 for the column order of unks
   (AnsatzSystem over one number field; the Expr route when a constant is not an
   algebraic number); {sub or None, neq} *)
SolveKernelSystem[T_, rem_, denv_, units_, unkLogs_, css_, monos_, gammas_, betas_, unks_, extra_] := Module[
  {sysK, q = T["q"], gens = T["gens"], nc = Length[css], V, EE, eqs, xs},
  sysK = Catch[AnsatzSystem[T, rem, denv, units, unkLogs, css, monos, gammas, betas, unks, extra], "PIM"];
  If[AssociationQ[sysK], Return[{sysK["sub"], sysK["neq"]}]];
  (* the Expr route: assembled without Cancel, row-reduced in the column order of unks *)
  V = PadRight[Table[Total[MapThread[#1 Times @@ (gens^#2) &, {css[[i]], monos}]]/denv, {i, nc}], T["n"]];
  EE = Together /@ (TowerD[T, V] - TPad[T, rem]);
  Do[EE = Together /@ (EE + gammas[[i]] TPad[T, AnalysisColumn[T, "u", units[[i, 1]]]]), {i, Length[units]}];
  Do[EE = Together /@ (EE + betas[[i]] TPad[T, AnalysisColumn[T, "s", unkLogs[[i, 1]]]]), {i, Length[unkLogs]}];
  Do[EE = Together /@ (EE + ex[[1]] TPad[T, ex[[2]]]), {ex, extra}];
  eqs = Flatten[If[# === 0, {}, CoefficientRules[Numerator[Together[#]], gens][[All, 2]]] & /@ EE];
  xs = LinSolveZero[eqs, unks];
  {If[xs === None, None, Thread[unks -> xs]], Length[eqs]}];

(* SolveWithKernels[key, split, retry, kernels, verbose, fixed]: Steps 15--20 of
   Part II's Algorithm 4 with the kernel columns (Algorithm S6, Steps 1--3):
   <|"I", "coeffs" (aligned with kernels), "bx" -> {bounds, exps, proved}|>, or
   {"nosol", bounds, proved} *)
SolveWithKernels[key_, split_, retry_, kernels_, verbose_, fixed_: None] := Module[{A = $analyses[key], T, f, Y, gens, q, m, nc,
    detLogs, denv, unitsBase, rem, rootLogs, unkLogs, sunits, bounds, exps, proved, bx, monos, css, units, gammas, betas, kappas,
    unks, extra, sub, neq, V, rat, y, surf, lrTerms, L, I0, coeffs},
  T = A["T"]; f = A["f"]; Y = A["Y"]; {gens, q, m} = {T["gens"], T["q"], T["m"]};
  nc = NCo[T];
  {detLogs, denv, unitsBase, rem, rootLogs} = Lookup[A, {"detLogs", "denv", "unitsBase", "rem", "rootLogs"}];
  {unkLogs, sunits} = AnalysisSpecials[key, split, verbose];
  {bounds, exps, proved} = If[fixed =!= None, fixed, ProxyBounds[key, unkLogs, kernels, retry, verbose]];
  bx = {bounds, exps, proved};
  Do[denv *= unkLogs[[i, 1]]^exps[[i]], {i, Length[unkLogs]}];
  monos = Tuples[Range[0, #] & /@ bounds];
  css = Table[cc[i - 1] @@@ monos, {i, nc}];
  units = Join[sunits, unitsBase];
  gammas = Table[gamma[i], {i, Length[units]}];
  betas = Table[beta[i], {i, Length[unkLogs]}];
  kappas = Table[kappa[i], {i, Length[kernels]}];
  unks = Join[Flatten[css], gammas, betas, kappas];
  extra = Table[{kappas[[j]], TPad[T, kernels[[j]]["col"]]}, {j, Length[kernels]}];
  {sub, neq} = SolveKernelSystem[T, rem, denv, units, unkLogs, css, monos, gammas, betas, unks, extra];
  If[verbose, Print["  ansatz: bounds ", bounds, ", ", Length[unks], " unknowns (", Length[kappas], " kernel columns)"]];
  If[sub === None, Return[{"nosol", bounds, proved}]];
  V = PadRight[Table[Total[MapThread[#1 Times @@ (gens^#2) &, {css[[i]], monos}]]/denv, {i, nc}], T["n"]];
  rat = Can /@ (V /. sub);
  y = If[q =!= None, q^(1/m), None];
  surf[u_] := If[q === None, u[[1]], Sum[TPad[T, u][[i + 1]] y^i/T["E"][[i + 1]], {i, 0, T["n"] - 1}]];
  lrTerms = Join[{#[[1]], ToY[T, #[[2]], Y]} & /@ detLogs,
    Table[{gammas[[i]] /. sub, ToY[T, units[[i, 1]], Y]}, {i, Length[units]}],
    Table[{betas[[i]] /. sub, unkLogs[[i, 1]]}, {i, Length[unkLogs]}]];
  {L, rat} = LogToReal[lrTerms, rat, f, T, Y, verbose];
  I0 = surf[rat] + If[q === None, L, L /. Y -> y] + Total[rootLogs];
  coeffs = kappas /. sub;
  <|"I" -> I0, "coeffs" -> coeffs, "bx" -> bx|>];

(* ------------------------------------------------------------- the driver (S7) *)
(* Part2[f, T, extended]: Part II's pipeline on f over T, with the published
   criteria or with (K5)-(K8), (B); no cached analysis is shared with another run *)
Part2[f_, T_, extended_: False] := With[{k = Hash[{CanT[f], T, TrueQ[extended]}]},
  If[AssociationQ[$part2Memo] && KeyExistsQ[$part2Memo, k], $part2Memo[k],
    With[{r = Block[{$analyses = <||>}, ExtendedBounds[extended, ParallelIntegrateMixed[f, T, "Verbose" -> False]]]},
      If[AssociationQ[$part2Memo], $part2Memo[k] = r]; r]]];
(* $part2Memo: Part II's results within one call of an entry point (the pipeline runs
   Part II on the same integrand more than once: E1, the certificates, E4); Part II is
   deterministic, so this saves time only *)
$part2Memo = None;
SetAttributes[WithPart2Memo, HoldAll];
WithPart2Memo[body_] := If[AssociationQ[$part2Memo], body, Block[{$part2Memo = <||>}, body]];

(* Certified[base]: Part II's status is a certificate of non-elementarity *)
Certified[base_] := ListQ[base] && Length[base] > 0 && base[[1]] === "not elementary";
NotTorsionQ[base_] := Certified[base] && Length[base] >= 2 && StringQ[base[[2]]] && StringContainsQ[base[[2]], "not torsion"];

Options[ParallelIntegrateSpecial] = {"Sources" -> None, "Back" -> None, "Y" -> None, "Verbose" -> False,
  "ElementaryFirst" -> True, "Strict" -> Automatic};

(* ParallelIntegrateSpecial[f, T, opts]: Algorithm S7.  Special functions are
   never introduced unnecessarily (Section 8.1): E1 an answer of Part II is
   returned unchanged; E2 Part II is rerun with the extended bounds before any
   kernel is introduced; E3 a third-kind kernel only for a residue divisor that
   Part II certifies non-torsion, an Ei kernel only for a non-constant residue;
   E4 every kernel removable at the same bounds is removed (greedily), and every
   special term that Part II integrates elementarily is replaced.  Necessity
   (Section 8.2): .certificate (Theorem 8.6, for f or the special part, over the
   tower or its pruned tower) and .independence (Theorem 8.7).  With strict a
   special answer only with a certificate, else {"failed", "special answer
   found, but the integrand is not certified non-elementary", answer}. *)
ParallelIntegrateSpecial[f00_List, T_Association, opts : OptionsPattern[]] := WithPart2Memo[Module[{r, f, cert, strict, verbose = OptionValue["Verbose"]},
  strict = If[OptionValue["Strict"] === Automatic, $StrictSP, TrueQ[OptionValue["Strict"]]];
  r = Pis[f00, T, OptionValue["Sources"], OptionValue["Back"], OptionValue["Y"], verbose, OptionValue["ElementaryFirst"]];
  If[! SpecialResultQ[r], Return[r]];
  f = TPad[T, f00];
  cert = If[r["certified"], r["certificate"], None];
  If[cert === None, cert = PMStage["CertifyNonelementary", CertifyNonelementary[f, T]]];
  If[cert === None,
    cert = PMStage["CertifyNonelementary", CertifyNonelementary[SpecialPart[r, T], T]];
    If[cert =!= None, cert = Prepend[cert, "special part"]]];
  r["certificate"] = cert; r["certified"] = cert =!= None;
  If[r["certified"],
    r["independence"] = If[Length[r["terms"]] == 1, {True, {"a single special term: necessary by the certificate"}},
      PMStage["IndependenceCertificate", IndependenceCertificate[r, T, r["base"], verbose]]];
    r["necessity"] = "necessary: the integrand has no elementary integral" <>
      If[r["independence"][[1]], "; the special terms are independent modulo elementary functions", ""],
    r["independence"] = {False, "no certificate for the integrand"};
    r["necessity"] = "minimal within the ansatz; non-elementarity not certified"];
  If[verbose, Print["  necessity: ", r["necessity"]]];
  If[strict && ! r["certified"], Return[{"failed", "special answer found, but the integrand is not certified non-elementary", r}, Module]];
  r]];

Pis[f00_, T_, sources0_, back_, Y_, verbose_, elementaryFirst_] := Catch[Module[{t0 = AbsoluteTime[], f, base = None, b2, certified,
    sources = sources0, determined = {}, f1, ei, cert, th, kernels, key, last = None, r, I0, used, mn},
  f = CanT[TPad[T, f00]];
  If[elementaryFirst,
    base = PMStage["Part2 (E1)", Part2[f, T]];
    If[ExprQ[base], If[verbose, Print["  elementary (Part II)"]]; Throw[base, "pis"]];     (* E1 *)
    If[verbose, Print["  Part II: ", StringTake[ToStr[base], UpTo[160]]]];
    If[! Certified[base],                                                                (* E2 *)
      b2 = PMStage["Part2 (E2 extended)", Part2[f, T, True]];
      If[ExprQ[b2], If[verbose, Print["  elementary (Part II with the extended bounds of Section 5)"]]; Throw[b2, "pis"]];
      If[Certified[b2], base = b2]]];
  certified = Certified[base];
  If[sources === None, sources = PMStage["ExpSources", ExpSources[T, back, Y, verbose]]];
  (* Step 2: Ei kernels determined by residues *)
  {f1, ei, cert} = PMStage["EiFromResidues", EiFromResidues[f, T, sources, verbose]];
  If[cert =!= None, Throw[cert, "pis"]];
  determined = Join[determined, ei];
  (* Step 3: third kind on a pencil, only for a certified non-torsion divisor   (E3) *)
  If[NotTorsionQ[base],
    {f1, th} = PMStage["ThirdKind", ThirdKind[f1, T, verbose]];
    determined = Join[determined, th]];
  (* Step 4: undetermined columns *)
  kernels = Join[PMStage["GammaCandidates", GammaCandidates[f1, T, sources, verbose]],
                 PMStage["EllipticColumns", EllipticColumns[T, verbose]]];
  (* Step 5: Part II's analysis of the remainder, the extended system *)
  If[AllZeroQ[f1], Throw[PMStage["Finish", Finish[SpecialResult[0, determined], T, verbose, t0, certified, base]], "pis"]];
  key = PMStage["SPAnalyse", SPAnalyse[f1, T, verbose]];
  If[ListQ[key],
    If[key[[1]] === "not elementary" && Length[key] == 3 && ! StringQ[key[[2]]],
      Throw[{"not in class", key[[2]], key[[3]], "non-constant residue of the remainder"}, "pis"]];
    Throw[{"failed", "analysis of the remainder", key}, "pis"]];
  Do[Do[
      r = PMStage["SolveWithKernels", SolveWithKernels[key, split, retry, kernels, verbose]];
      If[AssociationQ[r],
        mn = PMStage["Minimal", Minimal[key, split, retry, kernels, r["I"], r["coeffs"], r["bx"], verbose]];      (* E4 *)
        {I0, used} = mn;
        Throw[PMStage["Finish", Finish[SpecialResult[I0, Join[determined, used]], T, verbose, t0, certified, base]], "pis"]];
      last = r;
      If[r[[3]] && retry >= 0, Break[]],                  (* bounds proved: raising the guess changes nothing *)
      {retry, 0, 2}],
    {split, {False, True}}];
  {"failed", "no solution with the special kernels" <> If[last[[3]], " (every bound in force proved)", " (a guessed bound is in play)"],
    <|"bounds" -> last[[2]], "proved" -> last[[3]], "kernels" -> (KernelRepr /@ kernels), "part II" -> base|>}], "pis"];

(* Minimal: Algorithm S6, Steps 4--6 (E4): drop kernel columns greedily while the
   system, at the same bounds, stays solvable.  {I, {{coefficient, K}, ...}} *)
Minimal[key_, split_, retry_, kernels_, I0_, coeffs_, bx_, verbose_] := Module[{active, sol, changed = True, trial, r, cf},
  active = Select[Range[Length[kernels]], ! ExactZeroQ[coeffs[[#]]] &];
  sol = {I0, AssociationThread[active -> coeffs[[active]]]};
  While[changed && active =!= {},
    changed = False;
    Do[trial = DeleteCases[active, j];
      r = SolveWithKernels[key, split, retry, kernels[[trial]], False, bx];
      If[AssociationQ[r],
        cf = r["coeffs"];
        active = trial[[Select[Range[Length[trial]], ! ExactZeroQ[cf[[#]]] &]]];
        sol = {r["I"], AssociationThread[active -> cf[[Flatten[Position[trial, #] & /@ active]]]]};
        changed = True;
        If[verbose, Print["  kernel ", KernelRepr[kernels[[j]]], " removed: the system is solvable without it"]];
        Break[]],
      {j, active}]];
  {sol[[1]], Table[{sol[[2]][j], kernels[[j]]}, {j, active}]}];

(* SingleKernelNonelementaryQ[K, gens]: the column of this ONE kernel provably
   has no elementary antiderivative, so Part II cannot replace it and E4's probe
   on it is a search that cannot succeed.

   By construction -- the kernel identity Finish verifies just below, plus
   KernelAntiderivative -- a gamma kernel's column integrates to a nonzero
   constant multiple of

     s = 0                 ExpIntegralEi[v],
     s = a/k, k >= 2       Gamma[s, -v], or Sqrt[Pi] Erfc[h] when -v = h^2,

   and GammaCandidates builds s only as a/k with k >= 2 and GCD[a, k] == 1
   (:975), so s is NEVER a positive integer -- the one case in which Gamma[s, .]
   degenerates to an elementary polynomial times an exponential.  For a
   non-constant argument each of these is non-elementary over every elementary
   extension (Liouville, in Rosenlicht's form), which is the same theorem the
   stage's own certificates rest on.  Hence a single gamma term is never
   replaceable, and ONLY a cancellation between two or more kernels can make a
   sum of them elementary -- which is what the whole-part probe below tests, and
   the shape the coupling Lemma of AdditiveBlocks describes.

   Measured: E4 ran 344 full Part II integrations over the 312-case corpus and
   replaced nothing, 0 of 344, which is 14% of the stage's entire clock.  This
   is a soundness argument and not an empirical one: v non-constant is required
   explicitly, and the elliptic kinds are NOT claimed -- a degenerate pencil is
   ruled out upstream by BuildTower's squarefree radicand rather than here, so
   their probes are left in place. *)
SingleKernelNonelementaryQ[K_, gens_] := K["kind"] === "gamma" &&
  (K["s"] === 0 || ! IntegerQ[K["s"]]) && ! FreeQ[K["v"], Alternatives @@ gens];

(* ElementaryParts: Algorithm S7, Step 15 (E4): replace by Part II's elementary
   antiderivative the special part as a whole, else every special term that has one *)
ElementaryParts[res_, T_, verbose_] := Module[{nc = NCo[T], gens = T["gens"], col, whole, r, keep = {}, extra = 0},
  If[res["terms"] === {}, Return[res]];
  col[terms_] := Module[{u = TZero[T]}, Do[u = Padd[u, Pscale[t[[1]], TPad[T, t[[2]]["col"]]]], {t, terms}]; CanT[u]];
  whole = col[res["terms"]];
  If[AllTrue[Take[whole, nc], # === 0 &], Return[SpecialResult[res["elementary"], {}]]];
  (* the special part as a whole: with one term this is the per-term probe, so it
     is skipped on exactly the terms the invariant above covers *)
  If[Length[res["terms"]] >= 2 || ! SingleKernelNonelementaryQ[res["terms"][[1, 2]], gens],
    r = PMStage["Part2 (E4 whole)", Part2[whole, T, True]];
    If[ExprQ[r],
      PMStage["** E4 replaced the whole special part", Null];
      If[verbose, Print["  the special part integrates elementarily: replaced"]];
      Return[SpecialResult[res["elementary"] + r, {}]]]];
  Do[If[SingleKernelNonelementaryQ[t[[2]], gens], AppendTo[keep, t]; Continue[]];
    r = PMStage["Part2 (E4 per term)", Part2[col[{t}], T, True]];
    If[ExprQ[r],
      PMStage["** E4 replaced one term", Null];
      extra += r;
      If[verbose, Print["  ", KernelRepr[t[[2]]], " integrates elementarily: replaced"]],
      AppendTo[keep, t]],
    {t, res["terms"]}];
  SpecialResult[res["elementary"] + extra, keep]];

(* Finish: exact verification of every kernel identity, E4 (elementary parts
   replaced), and the necessity record *)
Finish[res0_, T_, verbose_, t0_, certified_: False, base_: None] := Module[{res = res0, K, vt, lhs, rhs},
  Do[K = t[[2]];
    If[K["kind"] === "gamma" && K["s"] =!= 0,
      vt = TScalar[T, K["v"]];
      lhs = TPad[T, TDiv[T, TowerD[T, K["omega"]], K["omega"]]];
      rhs = TPad[T, Padd[TowerD[T, vt], Pscale[K["s"] - 1, TDiv[T, TowerD[T, vt], vt]]]];
      If[! AllTrue[lhs - rhs, Can[#] === 0 &], Throw[{"exception", "kernel identity"}, "pis"]]],
    {t, res["terms"]}];
  res = ElementaryParts[res, T, verbose];
  If[res["terms"] === {},
    If[verbose, Print["  no special term is needed"]];
    Return[res["elementary"]]];
  res["certified"] = certified; res["base"] = base;
  res["certificate"] = If[certified, base, None];
  If[verbose, Print["  special answer: ", Length[res["terms"]], " special terms, ", ToString[NumberForm[AbsoluteTime[] - t0, {Infinity, 1}]], "s"]];
  res];

(* -------------------------------------------------------------- surface form (S9) *)
Num[e_, x_, x0_] := N[e /. x -> x0, 30];

(* the symbol of the radical among the back rules (the first Y-symbol, as the
   first Dummy of build_tower's back list) *)
BackY[back_, gens_] := SelectFirst[back[[All, 1]], MatchQ[#, _Symbol] && ! MemberQ[gens, #] &&
    StringMatchQ[SymbolName[#], "Y" ~~ ___] &, None];

SubstituteBack[e_, back_] := e //. back;

Options[IntegrateSurfaceSpecial] = {"Verbose" -> False, "Tower" -> None, "Samples" -> None,
  "Strict" -> None, "Details" -> False, "Additive" -> True};

(* IntegrateSurfaceSpecial[f, x, opts]: build the tower (BuildTower with the
   structure theorem, or "Tower" -> {T, fpair, back}), integrate with special
   functions, substitute back, fix the branch constant of every special term
   against its kernel numerically, and verify the whole answer by
   differentiation.  {answer, verified} or a status list. *)
IntegrateSurfaceSpecial[integrand_, x_Symbol, opts : OptionsPattern[]] := TimeConstrained[Module[{verbose = TrueQ[OptionValue["Verbose"]],
    tower = OptionValue["Tower"], integrandN = integrand /. PMExp -> Exp, bt, T, fp, back, Ysym, strict, res, samples, x0, yExpr, out, total, anti, kern, fixed, ver, blocks, sp},
  strict = If[OptionValue["Strict"] === None, $StrictSP, TrueQ[OptionValue["Strict"]]];
  (* the additive fast path.  Skipped when the caller pinned a tower (a tower
     for the whole integrand and a split of it are contradictory) and when
     "Details" is asked for, since the caller then wants this call's own
     internal result Association, which does not compose across blocks. *)
  If[tower === None && TrueQ[OptionValue["Additive"]] && ! TrueQ[OptionValue["Details"]],
    blocks = Quiet[AdditiveBlocks[integrand, x]];
    If[ListQ[blocks] && Length[blocks] >= 2,
      sp = SplitIntegrateSpecial[blocks, integrandN, x, strict, verbose, OptionValue["Samples"]];
      If[sp =!= None, Return[sp]]]];
  If[tower === None,
    bt = PMStage["BuildTower", Catch[Quiet[BuildTower[integrand, x, "StructureTheorem" -> True, "Verbose" -> verbose]], "build"]];
    If[! ListQ[bt] || Length[bt] != 4, Return[{"failed", "tower construction failed", ToStr[integrand]}]];
    {T, fp, back} = bt[[1 ;; 3]];
    If[verbose,
      Print["  tower: generators ", T["gens"], " with D = ", T["derivs"], If[T["q"] =!= None, ", y^" <> ToString[T["m"]] <> " = " <> ToStr[T["q"]], ""]];
      Print["  integrand: ", fp]],
    {T, fp, back} = tower];
  Ysym = BackY[back, T["gens"]];      (* strict was resolved above, for the split *)
  res = ParallelIntegrateSpecial[fp, T, "Back" -> back, "Y" -> Ysym, "Verbose" -> verbose, "Strict" -> False];
  samples = OptionValue["Samples"];
  If[samples === None,
    x0 = SamplePoint[integrandN, x];
    samples = {If[x0 =!= None, x0, 7/5], 9/4, 13/5}];
  PMStage["PresentSpecial", PresentSpecial[res, integrandN, x, T, back, samples, strict, verbose, TrueQ[OptionValue["Details"]]]]],
  $SpecialTimeBudget, {"failed", "time budget exceeded"}];

(* PresentSpecial: the surface form of a result of ParallelIntegrateSpecial: substitute
   back, fix the branch of every special term against its kernel, present, verify;
   {answer, verified} or the status (strict: an uncertified special answer is a status) *)
PresentSpecial[res_, integrandN_, x_, T_, back_, samples_, strict_, verbose_, details_] := Module[{yExpr, out, total, anti, kern, fixed, ver},
  yExpr = If[T["q"] =!= None, T["q"]^(1/T["m"]), None];
  If[ExprQ[res],
    out = SubstituteBack[res, back];
    Return[ConjugateIfBetter[out, integrandN, x, samples]]];
  If[! SpecialResultQ[res],
    Return[If[ListQ[res], Map[If[ExprQ[#] && ! StringQ[#], SubstituteBack[#, back], #] &, res], res]]];
  total = SubstituteBack[res["elementary"], back];
  Do[anti = SubstituteBack[KernelAntiderivative[t[[2]]], back];
    kern = t[[2]]["col"][[1]] + If[T["q"] =!= None, t[[2]]["col"][[2]] yExpr, 0];
    kern = SubstituteBack[kern, back];
    fixed = FixBranch[t[[1]] anti, t[[1]] kern, x, samples, t[[2]]];
    If[fixed === None,
      If[verbose, Print["  WARNING: branch of ", KernelRepr[t[[2]]], " not fixed numerically"]];
      fixed = t[[1]] anti];
    total += fixed,
    {t, res["terms"]}];
  total = HalfAngleFold[Present[total], back];
  (* the fold is INSIDE the verify: what the caller gets back is the form that
     was checked against the integrand, never a prettier unchecked cousin *)
  {total, ver} = ConjugateIfBetter[total, integrandN, x, samples];
  If[strict && ! res["certified"], Return[{"failed", "special answer found, but the integrand is not certified non-elementary", total}]];
  If[details, {total, ver, res}, {total, ver}]];

(* FixBranch[anti, kern, x, samples, K]: anti up to a root of unity / sign: the
   multiplier making d anti / dx = kern at the sample points *)
FixBranch[anti_, kern_, x_, samples_, K_] := Module[{mults = {1, -1, I, -I}, k, d, ok, a, b},
  If[K["kind"] === "gamma" && K["s"] =!= 0,
    k = Denominator[K["s"]];
    mults = Join[mults, Table[Exp[2 Pi I j/k], {j, 1, k - 1}]]];
  d = QuietCheck[D[anti, x], $Failed];
  If[d === $Failed, Return[None]];
  Do[ok = True;
    Do[a = Quiet[Num[mu d, x, x0]]; b = Quiet[Num[kern, x, x0]];
      If[! NumericQ[a] || ! NumericQ[b] || Abs[a - b] > 10^-15 (1 + Abs[b]), ok = False; Break[]],
      {x0, samples}];
    If[ok, Return[mu anti, Module]],
    {mu, mults}];
  None];

(* ConjugateIfBetter[e, integrand, x, samples]: {answer, verified}, trying the
   OTHER embedding of the constant field when the first does not verify.  The
   ansatz is solved over Q(i) (or a field containing it) and nothing in the
   linear algebra chooses between i -> +-i: both solve, and only one
   differentiates back once the surface radicals are read on their principal
   branches.  Conjugating the coefficients leaves the functions alone, so this is
   not reachable by the kernel-wise multiplier FixBranch applies.  Sound: the
   conjugate is one more CANDIDATE, accepted only on the same evidence as the
   original -- VerifyAnswer at every sample.  Corpus #34,
   ArcSin[Sqrt[x+1]]/Sqrt[x], whose surface is 2 Sqrt[x] ArcSin[Sqrt[1+x]]
   - 2 I Sqrt[1+x] where the antiderivative has +2 I Sqrt[1+x]. *)
ConjugateIfBetter[e_, integrand_, x_, samples_] := Module[{ver, ec},
  ver = VerifyAnswer[e, integrand, x, samples];
  If[ver || FreeQ[e, _Complex], Return[{e, ver}]];
  ec = e /. Complex[re_, im_] :> Complex[re, -im];
  If[TrueQ[VerifyAnswer[ec, integrand, x, samples]], {ec, True}, {e, ver}]];

VerifyAnswer[expr_, integrand_, x_, samples_] := Module[{d, a, b},
  d = QuietCheck[D[expr, x], $Failed];
  If[d === $Failed, Return[False]];
  AllTrue[samples, (a = Quiet[Num[d, x, #]]; b = Quiet[Num[integrand, x, #]];
    NumericQ[a] && NumericQ[b] && Abs[a - b] < 10^-15 (1 + Abs[b])) &]];

(* ------------------------------------------ bounds at the places (Section 5, S5) *)
(* Part II decides the hypotheses of its Theorem 8.7 at a place by the criteria
   (K1)-(K4) of Proposition 8.11 and leaves every other place to the classical
   guess.  The extensions below, proved in Section 5 of the paper, decide more
   places; they are installed for the special stage only (ExtendedBounds):
     (K5) extra constants (Theorem 5.1, Proposition 5.2);
     (K6) a rational residue field, one attaining generator;
     (K7) two primitive layers;
     (K8) exponentials over primitives;
     (B)  branch special primes of a square-root radical (the uniformiser y).
   The names of the extended constants: {"rho", n} the class of D pi / pi^n,
   {"D", t_k, n} the class of D t_k / pi^n, {"rhoexpr", n}, {"Dexpr", t_k, n}
   the same as expressions in the other generators.                         *)

(* BranchClass[T, e, pp, g, n, gens]: the class of the scalar e / y^n at the
   ramified place over pp (m = 2, y^2 = q = pp q1), when it is a constant *)
BranchClass[T_, e_, pp_, g_, n_, gens_] := If[EvenQ[n], ClassMod[Can[e/Can[T["q"]/pp]^(n/2)], pp, g, n/2, gens], None];

(* ClassExpr[e, pp, g, n, gens]: the class of e / pp^n modulo pp as an expression
   in the other generators (pp linear in g), or None *)
ClassExpr[e_, pp_, g_, n_, gens_] := QuietCheck[Module[{cl},
    If[Vp[e, pp, gens] < n, Return[None, Module]];
    If[! PolynomialQ[pp, g] || Exponent[pp, g] != 1, Return[None, Module]];
    cl = CoefficientList[pp, g];
    Can[Can[e/pp^n] /. g -> -cl[[1]]/cl[[2]]]], None];

SpecialConstsExt[T_, pp_, Dpi_, s_, g_, names_] := Module[{gens = T["gens"], q = T["q"], branch, out = <||>, kind, n, u, a, b, q1, cls, ok},
  branch = q =!= None && T["m"] == 2 && Vp[q, pp, gens] > 0;
  Do[
    If[ListQ[nm],
      kind = nm[[1]]; n = Last[nm];
      u = TPad[T, If[MemberQ[{"rho", "rhoexpr"}, kind], Dpi, T["derivs"][[Position[gens, nm[[2]]][[1, 1]]]]]];
      If[kind === "Dexpr",
        out[nm] = Which[
          ! branch && AllTrue[Rest[u], Can[#] === 0 &], ClassExpr[u[[1]], pp, g, n, gens],
          branch && EvenQ[n] && AllTrue[Rest[u], Can[#] === 0 &], ClassExpr[Can[u[[1]]/Can[q/pp]^(n/2)], pp, g, n/2, gens],
          True, None];
        Continue[]];
      If[kind === "rhoexpr",
        {a, b} = u[[1 ;; 2]];
        out[nm] = Which[
          ! branch, If[AllTrue[Rest[u], Can[#] === 0 &], ClassExpr[a, pp, g, n, gens], None],
          OddQ[n] && (Can[a] === 0 || 2 Vp[a, pp, gens] > n), q1 = Can[q/pp]; ClassExpr[Can[b/q1^((n - 1)/2)], pp, g, (n - 1)/2, gens],
          True, None];
        Continue[]],
      If[nm === "pi", u = TPad[T, Dpi]; n = s + 1,
        u = Can[#/nm] & /@ TPad[T, T["derivs"][[Position[gens, nm][[1, 1]]]]]; n = s]];
    If[! branch,
      cls = If[Can[#] === 0, 0, ClassMod[#, pp, g, n, gens]] & /@ u;
      out[nm] = If[cls[[1]] =!= None && AllTrue[Rest[cls], # === 0 &], cls[[1]], None];
      Continue[]];
    (* branch: pi = y; an element a + b y has valuation min(2 v(a), 2 v(b) + 1) *)
    {a, b} = u[[1 ;; 2]];
    If[EvenQ[n],
      ok = Can[b] === 0 || 2 Vp[b, pp, gens] + 1 > n;
      out[nm] = If[ok && Can[a] =!= 0, BranchClass[T, a, pp, g, n, gens], If[ok, 0, None]],
      ok = Can[a] === 0 || 2 Vp[a, pp, gens] > n;
      out[nm] = If[ok && Can[b] =!= 0, BranchClass[T, b, pp, g, n - 1, gens], If[ok, 0, None]]],
    {nm, names}];
  out];

InfConstsExt[T_, g_, pl_, names_] := Module[{gens = T["gens"], plain, out, u, lc},
  plain = Select[names, ! ListQ[#] &];
  out = If[plain =!= {}, Association[Normal[InfConstsOrig[T, g, pl, plain]]], <||>];
  Do[If[! ListQ[nm], Continue[]];
    u = If[MemberQ[{"rho", "rhoexpr"}, nm[[1]]], -T["derivs"][[Position[gens, g][[1, 1]]]], T["derivs"][[Position[gens, nm[[2]]][[1, 1]]]]];
    lc = LcPlace[T, u, g];
    out[nm] = If[lc =!= None && pl <= Length[lc], lc[[pl]], None],
    {nm, names}];
  out];

ConstValQ[v_] := NumericQ[v];
(* LexLess[a, b]: the tuple a precedes b lexicographically (Python's tuple order) *)
LexLess[a_, b_] := Module[{k = 1}, While[k <= Min[Length[a], Length[b]] && a[[k]] == b[[k]], k++];
  If[k > Min[Length[a], Length[b]], Length[a] < Length[b], a[[k]] < b[[k]]]];

(* DecideBoundExt: Part II's criteria, then (K5)-(K8) of Proposition 5.2 *)
DecideBoundExt[T_, own_, mono_, r_, A_, kinds_, rhoFree_, lamK1_, lamK4_, rhoR_, upperAllAttain_, consts_, ext_: None] := Module[
  {c, s, shifts, sigPi, lam, gens = T["gens"], rr, AA, BB, names, vals, rho, mus},
  c = DecideBoundOrig[T, own, mono, r, A, kinds, rhoFree, lamK1, lamK4, rhoR, upperAllAttain, consts, ext];
  If[c =!= None || consts === None || ext === None, Return[c]];
  {s, shifts, sigPi, lam} = Lookup[ext, {"s", "shifts", "sigPi", "lam"}];
  rr = sigPi - s;
  If[rr == 0, Return[DecideK6[T, own, s, shifts, lam, consts]]];
  If[KeyExistsQ[shifts, "y"] && shifts["y"] <= s + rr, Return[None]];
  AA = Select[gens, Lookup[shifts, #, None] === s && ! MemberQ[own, #] &];
  BB = Select[gens, ! MemberQ[own, #] && ! MemberQ[AA, #] &];
  If[AnyTrue[AA, Kind[T, Position[gens, #][[1, 1]]] =!= "exp" &], Return[None]];
  If[AnyTrue[AA, Intersection[lam[#], AA] =!= {} &], Return[None]];
  If[AnyTrue[BB, Lookup[shifts, #, Infinity] < s + rr &], Return[None]];
  names = Join[{{"rho", s + 1 + rr}}, AA, {"D", #, s + rr} & /@ Select[BB, Lookup[shifts, #, Infinity] === s + rr &]];
  vals = QuietCheck[consts[names], $Failed];
  If[! AssociationQ[vals] || AnyTrue[Values[vals], # === None &], Return[None]];
  If[! AllTrue[Values[vals], ConstValQ], Return[None]];
  rho = vals[{"rho", s + 1 + rr}];
  If[ExactZeroQ[rho], Return[None]];
  mus = vals[#] & /@ AA;
  If[AnyTrue[mus, ExactZeroQ], Return[None]];
  Do[If[InQSpanQ[mus[[i]], Delete[mus, i]], Return[None, Module]], {i, Length[mus]}];
  "K5"];

(* DecideK6: (K6), (K7), (K8) of Proposition 5.2, all with r = 0 *)
DecideK6[T_, own_, s_, shifts_, lam_, consts_] := Module[{gens = T["gens"], AA, c, vals, rho, lamc, R, n, d, rts, res, k8, l1, l2, R2},
  AA = Select[gens, Lookup[shifts, #, None] === s && ! MemberQ[own, #] &];
  If[KeyExistsQ[shifts, "y"] && shifts["y"] === s, Return[None]];
  If[Length[AA] == 1,
    c = AA[[1]];
    If[Kind[T, Position[gens, c][[1, 1]]] === "exp", Return[None]];
    vals = QuietCheck[consts[{{"rhoexpr", s + 1}, {"Dexpr", c, s}}], $Failed];
    If[! AssociationQ[vals] || AnyTrue[Values[vals], # === None &], Return[None]];
    {rho, lamc} = {vals[{"rhoexpr", s + 1}], vals[{"Dexpr", c, s}]};
    If[ExactZeroQ[lamc] || ! SubsetQ[{c}, FreeSyms[lamc]] || ! SubsetQ[{c}, FreeSyms[rho]], Return[None]];
    R = Can[rho/lamc];
    If[R === 0, Return[None]];
    {n, d} = {Numerator[R], Denominator[R]};
    If[Exponent[n, c] >= Exponent[d, c], Return["K6"]];                  (* a polynomial part: never a log-derivative *)
    If[AnyTrue[FactorList[d], ! FreeQ[#[[1]], c] && #[[2]] >= 2 &], Return["K6"]];   (* a multiple pole *)
    rts = QuietCheck[RootsOf[d, c], $Failed];
    If[rts === $Failed, Return[None]];
    res = Simplify[(n /. c -> #)/(D[d, c] /. c -> #)] & /@ rts;
    If[Length[res] != Exponent[d, c], Return[None]];
    If[AllTrue[res, MatchQ[#, _Integer | _Rational] &], Return[None]];  (* -a rho / lambda is a log-derivative for some a < 0 *)
    Return["K6"]];
  k8 = DecideK8[T, AA, s, consts];
  If[k8 =!= None, Return[k8]];
  If[Length[AA] == 2 && AllTrue[AA, Kind[T, Position[gens, #][[1, 1]]] === "prim" &],
    vals = QuietCheck[consts[Join[{{"rhoexpr", s + 1}}, {"Dexpr", #, s} & /@ AA]], $Failed];
    If[! AssociationQ[vals] || AnyTrue[Values[vals], # === None &], Return[None]];
    rho = vals[{"rhoexpr", s + 1}];
    If[ExactZeroQ[rho] || FreeSyms[rho] =!= {}, Return[None]];
    Do[{l1, l2} = {vals[{"Dexpr", pr[[1]], s}], vals[{"Dexpr", pr[[2]], s}]};
      If[ExactZeroQ[l1] || FreeSyms[l1] =!= {} || ! SubsetQ[{pr[[1]]}, FreeSyms[l2]] || ExactZeroQ[l2], Continue[]];
      R2 = Can[l2/l1];
      d = Denominator[R2];
      If[FreeQ[d, pr[[1]]], Continue[]];
      rts = QuietCheck[RootsOf[d, pr[[1]]], $Failed];
      If[rts === $Failed, Continue[]];
      If[AnyTrue[rts, ! ExactZeroQ[Quiet[Residue[R2, {pr[[1]], #}]]] &], Return["K7", Module]],
      {pr, {AA, Reverse[AA]}}]];
  None];

(* DecideK8: (K8), r = 0: the attaining generators are primitives with constant
   leading coefficients and hyperexponentials whose lambda is a polynomial in
   those primitives, and rho is a nonzero constant: 'K8' when no rational
   combination of the lambdas is a nonzero constant *)
DecideK8[T_, AA_, s_, consts_] := Module[{gens = T["gens"], prims, exps, vals, rho, lams = {}, lam, ns, comb, eqs, const, sol, cval, rules},
  If[AA === {}, Return[None]];
  prims = Select[AA, Kind[T, Position[gens, #][[1, 1]]] === "prim" &];
  exps = Select[AA, Kind[T, Position[gens, #][[1, 1]]] === "exp" &];
  If[Length[prims] + Length[exps] != Length[AA] || exps === {}, Return[None]];
  vals = QuietCheck[consts[Join[{{"rhoexpr", s + 1}}, {"Dexpr", #, s} & /@ AA]], $Failed];
  If[! AssociationQ[vals] || AnyTrue[Values[vals], # === None &], Return[None]];
  rho = vals[{"rhoexpr", s + 1}];
  If[ExactZeroQ[rho] || FreeSyms[rho] =!= {}, Return[None]];
  Do[If[FreeSyms[vals[{"Dexpr", g, s}]] =!= {} || ExactZeroQ[vals[{"Dexpr", g, s}]], Return[None, Module]], {g, prims}];
  Do[lam = Can[vals[{"Dexpr", g, s}]/g];
    If[If[prims =!= {}, ! SubsetQ[prims, FreeSyms[lam]] || ! PolynomialQ[lam, prims], FreeSyms[lam] =!= {}], Return[None, Module]];
    AppendTo[lams, Expand[lam]], {g, exps}];
  ns = Table[Unique["n"], {Length[lams]}];
  comb = Expand[ns . lams];
  If[prims =!= {},
    rules = CoefficientRules[comb, prims];
    eqs = Select[rules, AnyTrue[#[[1]], # != 0 &] &][[All, 2]];
    const = Total[Select[rules, AllTrue[#[[1]], # == 0 &] &][[All, 2]]],
    eqs = {}; const = comb];
  sol = If[eqs === {}, {{}}, Quiet[Solve[Thread[eqs == 0], ns]]];
  If[sol === {}, Return["K8"]];
  cval = Expand[const /. First[sol]];
  If[cval === 0, "K8", None]];

(* SpecialDataExt: Part II's SpecialData, extended to the ramified place over a
   special prime of a square-root radical (the uniformiser y, exact valuations) *)
SpecialDataExt[key_, T_, pp_] := Module[{gens = T["gens"], q = T["q"], m = T["m"], res, parents, par, ownRes, own, vp, shifts, vD,
    pi, Dpi, sigPi, s, AA, lam, lamSyms, rhoFree, consts, crit},
  If[! (q =!= None && m == 2 && Vp[q, pp, gens] > 0),
    res = SpecialDataOrig[key, T, pp];
    parents = $analyses[key]["parent"];
    If[res[[1]] === None && KeyExistsQ[parents, pp],
      (* a factor split off a special over Fbar: decide at the factor itself (its
         conjugates carry the same decision) rather than at the parent *)
      par = parents[pp];
      $analyses[key, "parent"] = KeyDrop[parents, {pp}];
      $analyses[key, "spec"] = KeyDrop[$analyses[key]["spec"], {pp}];
      ownRes = SpecialDataOrig[key, T, pp];
      $analyses[key, "parent", pp] = par;
      If[ownRes[[1]] =!= None, $analyses[key, "spec", pp] = ownRes; Return[ownRes]];
      $analyses[key, "spec", pp] = res];                  (* no decision at the factor: keep the parent's record *)
    Return[res]];
  If[KeyExistsQ[$analyses[key]["spec"], pp], Return[$analyses[key]["spec"][pp]]];
  own = Select[gens, ! FreeQ[pp, #] &];
  vp[u_] := VP[u, pp, T, True];
  shifts = <||>;
  Do[vD = vp[T["derivs"][[k]]];
    If[vD =!= Infinity, shifts[gens[[k]]] = vD - vp[TScalar[T, gens[[k]]]]], {k, Length[gens]}];
  pi = TUnit[T, 1];
  Dpi = TowerD[T, pi];
  sigPi = vp[Dpi] - vp[pi];
  s = Min[Append[Values[shifts], sigPi]];
  AA = Select[gens, Lookup[shifts, #, None] === s && ! MemberQ[own, #] &];
  lam = Table[LamCoeffs[T, k, None], {k, Length[gens]}];
  lamSyms = AssociationThread[gens -> Table[Complement[Select[gens, Function[z, ! FreeQ[lam[[k]], z]]], own], {k, Length[gens]}]];
  rhoFree = FreeQ[DeleteCases[Can /@ Dpi, 0], Alternatives @@ Append[AA, Unique[]]];
  consts = If[Length[own] == 1, Function[names, SpecialConstsExt[T, pp, Dpi, s, own[[1]], names]], None];
  crit = DecideBound[T, own, None, sigPi - s, AA, Kind[T, Position[gens, #][[1, 1]]] & /@ AA, rhoFree,
    AllTrue[AA, FreeQ[lam[[Position[gens, #][[1, 1]]]], Alternatives @@ Append[AA, Unique[]]] &],
    AllTrue[AA, FreeQ[lam[[Position[gens, #][[1, 1]]]], Alternatives @@ Append[Complement[gens, own], Unique[]]] &],
    None, AllTrue[Complement[gens, own], MemberQ[AA, #] &], consts,
    <|"s" -> s, "shifts" -> shifts, "sigPi" -> sigPi, "lam" -> lamSyms, "own" -> own|>];
  $analyses[key, "spec", pp] = {crit, s, sigPi - s, m, True};
  $analyses[key]["spec"][pp]];

(* ExtendedBounds[on, body]: (K5)-(K8) and the branch special primes (B) installed
   into Part II's Algorithm 6 for the evaluation of body (a Block over the four
   hooks: restored on exit, Throw or Abort) *)
SetAttributes[ExtendedBounds, HoldRest];
ExtendedBounds[on_, body_] := If[! TrueQ[on], body,
  Block[{DecideBound, SpecialData, InfConsts, SpecialConsts},
    DecideBound[args___] := DecideBoundExt[args];
    SpecialData[args___] := SpecialDataExt[args];
    InfConsts[args___] := InfConstsExt[args];
    SpecialConsts[args___] := SpecialConstsExt[args];
    body]];

(* ------------------------------------------ necessity certificates (Section 8.2, S8) *)
(* PrunedTower[T, f]: Lemma 8.5, the sub-tower L0 generated by the generators
   that f, q and the derivatives of the kept generators involve, when L is an
   elementary extension of L0 (every removed generator a logarithm c log u or an
   exponential e^v of L0, by Part II); a pruned tower with a single generator g,
   D g = r(g) != 1, is rescaled to d/dg and f divided by r.  {T0, f0} or None *)
PrunedTower[T_, f_] := Catch[Module[{gens = T["gens"], q = T["q"], need, changed = True, new, kept, T0, d, kind, r, lgs, Ld, cL, ok, f0, g0, d0, r0},
  need = Select[gens, ! FreeQ[f, #] &];
  If[q =!= None, need = Union[need, Select[gens, ! FreeQ[q, #] &]]];
  While[changed,
    changed = False;
    Do[If[MemberQ[need, gens[[k]]],
        new = Complement[Select[gens, Function[z, ! FreeQ[T["derivs"][[k]], z]]], need];
        If[new =!= {}, need = Union[need, new]; changed = True]],
      {k, Length[gens]}]];
  kept = Select[gens, MemberQ[need, #] &];
  If[Length[kept] == Length[gens] || kept === {}, Throw[None, "pt"]];
  T0 = Tower[kept, T["derivs"][[Position[gens, #][[1, 1]]]] & /@ kept, q, T["m"]];
  Do[If[MemberQ[need, gens[[k]]], Continue[]];
    d = T["derivs"][[k]];
    If[AnyTrue[Rest[d], Can[#] =!= 0 &] || Complement[FreeSyms[d[[1]]], need, {gens[[k]]}] =!= {}, Throw[None, "pt"]];
    kind = Kind[T, k];
    ok = QuietCheck[Which[
      kind === "prim",
        r = Block[{$analyses = <||>}, ParallelIntegrateMixed[TScalar[T0, d[[1]]], T0]];
        (* one logarithm c log u; an arctangent and an artanh are single logarithms as well *)
        lgs = DeleteDuplicates[Cases[{r}, _Log | _ArcTan | _ArcTanh | _RootSum, Infinity]];
        If[ExprQ[r] && Length[lgs] <= 1 && FreeQ[lgs, RootSum],
          If[lgs === {}, True,
            Ld = Unique["L"];
            cL = D[Expand[r] /. lgs[[1]] -> Ld, Ld];
            FreeQ[cL, Alternatives @@ Append[gens, Ld]]],
          False],
      kind === "exp",
        r = Block[{$analyses = <||>}, ParallelIntegrateMixed[TScalar[T0, Can[d[[1]]/gens[[k]]]], T0]];
        ExprQ[r] && FreeQ[r, _Log | _ArcTan | _ArcTanh | _RootSum],
      True, False], False];
    If[! TrueQ[ok], Throw[None, "pt"]],
    {k, Length[gens]}];
  f0 = TPad[T0, f];
  If[Length[kept] == 1 && T0["q"] =!= None,
    g0 = kept[[1]]; d0 = T0["derivs"][[1]]; r0 = d0[[1]];
    If[r0 =!= 1 && AllTrue[Rest[d0], Can[#] === 0 &] && SubsetQ[{g0}, FreeSyms[r0]],
      T0 = Tower[{g0}, {{1}}, T0["q"], T0["m"]];
      f0 = Can[#/r0] & /@ f0]];
  {T0, f0}], "pt"];

(* CertifyNonelementary[f, T]: a Part II certificate that f has no elementary
   integral, over T or over its pruned tower, with the published and with the
   extended criteria; the certificate or None *)
CertifyNonelementary[f_, T_] := Module[{pr = PrunedTower[T, f], r},
  Do[Do[r = QuietCheck[Part2[tf[[2]], tf[[1]], ext], $Failed];
      If[Certified[r], Return[If[tf[[1]] === T, r, Join[{"pruned tower", tf[[1]]["gens"]}, r]], Module]],
      {ext, {False, True}}],
    {tf, Join[{{T, f}}, If[pr =!= None, {pr}, {}]]}];
  None];

SpecialPart[res_, T_] := Module[{u = TZero[T]},
  Do[u = Padd[u, Pscale[t[[1]], TPad[T, t[[2]]["col"]]]], {t, res["terms"]}];
  CanT[u]];

(* IndependenceCertificate[res, T, base, verbose]: Algorithm S8, Steps 5--12
   (Theorem 8.7): no nontrivial constant combination of the kernels of res has an
   elementary integral over L.  {True, reasons} or {False, reason} *)
IndependenceCertificate[res_, T_, base_, verbose_: False] := Catch[Module[{Ks, places, iEi, iTh, iN, reasons = {}, taus, cs, eqs, e, b, M, ok, why, ok2, why2, others},
  Ks = res["terms"][[All, 2]];
  (* an Ei column with a normal zero carries a residue there *)
  places = Table[If[K["kind"] === "gamma" && K["s"] === 0, If[KeyExistsQ[K, "place"], K["place"], NormalZero[T, K]], None], {K, Ks}];
  iEi = Select[Range[Length[Ks]], Ks[[#]]["kind"] === "gamma" && Ks[[#]]["s"] === 0 && places[[#]] =!= None &];
  iTh = Select[Range[Length[Ks]], Ks[[#]]["kind"] === "third" &];
  iN = Complement[Range[Length[Ks]], iEi, iTh];
  (* (a) Ei kernels found by residues: 1 and their residues independent at every place *)
  Do[taus = Table[ResidueAt[T, Ks[[j]]["col"], pl[[1]], pl[[2]], pl[[1]] - pl[[2]]], {j, iEi}];
    taus = Select[taus, ! AllTrue[#, Can[#] === 0 &] &];
    cs = Table[Unique["c"], {Length[taus] + 1}];
    eqs = {};
    Do[e = If[i == 1, cs[[1]], 0] + Sum[cs[[j + 1]] taus[[j, i]], {j, Length[taus]}];
      eqs = Join[eqs, NumCoeffs[e, T["gens"]]], {i, NCo[T]}];
    ok = If[eqs === {}, False, {b, M} = Normal[CoefficientArrays[eqs, cs]]; MatrixRank[M] == Length[cs]];
    If[! ok, Throw[{False, "residues of the Ei kernels dependent at " <> ToStr[pl[[1]]] <> " = " <> ToStr[pl[[2]]]}, "ic"]];
    AppendTo[reasons, "Ei residues independent at " <> ToStr[pl[[1]]] <> " = " <> ToStr[pl[[2]]]],
    {pl, DeleteDuplicates[places[[iEi]]]}];
  (* (b) at most one third-kind kernel, with a divisor certified non-torsion *)
  If[Length[iTh] > 1, Throw[{False, "several third-kind kernels: the joint divisor is not certified"}, "ic"]];
  If[iTh =!= {},
    If[! NotTorsionQ[base], Throw[{False, "third-kind divisor not certified non-torsion"}, "ic"]];
    AppendTo[reasons, "third-kind divisor certified non-torsion (Part II, Prop. 9.4)"]];
  (* (c) the residue-free kernels: the parametric holomorphic-remainder certificate *)
  Do[others = Ks[[DeleteCases[iN, j]]];
    {ok, why} = ParametricCertificate[T, Ks[[j]], others];
    If[! ok,
      {ok2, why2} = ParametricPruned[T, Ks[[j]], others];
      If[! ok2, Throw[{False, KernelRepr[Ks[[j]]] <> ": " <> why}, "ic"]]],
    {j, iN}];
  If[iN =!= {}, AppendTo[reasons, ToString[Length[iN]] <> " residue-free kernels: parametric certificate"]];
  {True, reasons}], "ic"];

(* NormalZero[T, K]: a place {g, rho} over a normal prime at which the s = 0
   kernel K has a non-constant residue, or None *)
NormalZero[T_, K_] := Module[{num = Numerator[Together[K["v"]]], syms, cl, pls, tau},
  syms = Select[T["gens"], ! FreeQ[num, #] &];
  If[syms === {}, Return[None]];
  Do[cl = QuietCheck[ClassifyPrime[T, fac[[1]]], $Failed];
    If[cl === $Failed || cl[[4]] || cl[[1]], Continue[]];
    pls = LinearPlaces[T, fac[[1]]];
    Do[tau = ResidueAt[T, TPad[T, K["col"]], pl[[1]], pl[[2]], pl[[3]]];
      If[! IsConst[T, tau, pl[[1]], pl[[2]]], Return[pl[[1 ;; 2]], Module]],
      {pl, If[pls === None, {}, pls]}],
    {fac, FactorsOf[num]}];
  None];

(* ParametricPruned: the parametric certificate over the pruned (and rescaled)
   tower of Lemma 8.5, all kernels lying in L0 *)
ParametricPruned[T_, Kj_, others_] := Module[{allk = Prepend[others, Kj], pr, T0, kept, g0, r0, moved, ok, why},
  pr = PrunedTower[T, SpecialPart[SpecialResult[0, Table[{3 + 2 (i - 1), allk[[i]]}, {i, Length[allk]}]], T]];
  If[pr === None, Return[{False, "no pruned tower"}]];
  T0 = pr[[1]]; kept = T0["gens"];
  If[AnyTrue[allk, HasAny[#["col"], Complement[T["gens"], kept]] &], Return[{False, "kernels outside the pruned tower"}]];
  g0 = kept[[1]];
  r0 = If[Length[kept] == 1 && T0["derivs"][[1, 1]] === 1, T["derivs"][[Position[T["gens"], g0][[1, 1]], 1]], 1];
  moved = MakeKernel[#["kind"], Can[#/r0] & /@ TPad[T0, #["col"]]] & /@ allk;
  {ok, why} = ParametricCertificate[T0, moved[[1]], Rest[moved]];
  {ok, why <> " (pruned tower, Lemma 8.5)"}];

(* ParametricCertificate[T, Kj, others]: no combination K_j - Sum c_i K_i has an
   elementary integral: the system with right side K_j and the other kernels as
   columns has no solution at bounds that are all proved, and the other
   hypotheses of Part II's Proposition 9.2(b) hold for the family *)
ParametricCertificate[T_, Kj_, others_] := Module[{col = CanT[TPad[T, Kj["col"]]], key, A, gens = T["gens"], q = T["q"], unk0, splittable,
    split, unkLogs, curve, bounds, exps, proved, r, deg2},
  key = SPAnalyse[col, T, False];
  If[ListQ[key], Return[{False, "analysis: " <> StringTake[ToStr[key], UpTo[60]]}]];
  A = $analyses[key];
  deg2[pp_] := AnyTrue[gens, PolynomialQ[pp, #] && Exponent[pp, #] >= 2 &];
  unk0 = AnalysisSpecials[key, False, False][[1]];
  splittable = AnyTrue[unk0, deg2[#[[1]]] &];
  split = splittable;
  unkLogs = AnalysisSpecials[key, split, False][[1]];
  curve = Select[gens, q =!= None && ! FreeQ[q, #] &];
  If[AnyTrue[unkLogs, q =!= None && SubsetQ[curve, FreeSyms[#[[1]]]] &],
    Return[{False, "a special over the curve variable (S'-units searched within a bound)"}]];
  If[split && AnyTrue[unkLogs, deg2[#[[1]]] &], Return[{False, "a special irreducible of degree >= 2 over Fbar"}]];
  {bounds, exps, proved} = ProxyBounds[key, unkLogs, others, 0, False];
  If[! proved, Return[{False, "a guessed bound is in play"}]];
  If[! (A["unitsComplete"] || TrueQ[QuietCheck[SecondKindAtInfinityQ[T, A["rem"]], False]]), Return[{False, "unit candidates not complete"}]];
  r = SolveWithKernels[key, split, 0, others, False, {bounds, exps, proved}];
  If[AssociationQ[r], Return[{False, "dependent: the system has a solution"}]];
  If[! SPVerifiedResidueFree[A["rem"], T], Return[{False, "residual not verified residue-free"}]];
  Do[If[! SPVerifiedResidueFree[CanT[TPad[T, K["col"]]], T], Return[{False, "a kernel column not verified residue-free"}, Module]], {K, others}];
  {True, "parametric certificate"}];

SPVerifiedResidueFree[rem_, T_] := TrueQ[QuietCheck[
  If[Length[T["gens"]] == 1 && T["q"] =!= None && T["m"] == 2 && T["derivs"][[1]] === {1, 0},
    ResidueFree[T, rem, Unique["y"]],
    Block[{$analyses = <||>}, VerifiedResidueFree[rem, T]]], False]];

(* ------------------------------------------------ partial answers (Section 8.3, S10) *)
(* When no complete antiderivative is found, PartialIntegrate returns an exact
   decomposition f = D(I) + r with I elementary (or elementary plus certified
   special terms) and a REDUCED remainder r: the Hermite part is always
   integrated; a non-constant residue outside the exponential span is carried by
   r; a residue divisor certified non-torsion by the third-kind kernel, in r
   (elementary mode) or as Pi (special mode); among such remainders r is chosen
   by a reduction with the ansatz columns before the remainder columns (the
   additive terms of f, then the reduced monomial basis, simplest first). *)

(* ResidueSplits[f, T, sources, special, verbose]: Step (P1): at every normal
   prime with a critical (or deeper) pole and a non-constant residue tau: in
   special mode subtract the Ei kernels of the part of tau in the exponential
   span; move what is left of tau to the remainder.  {f', remainder, Ei terms} *)
ResidueSplits[f0_, T_, sources_, special_, verbose_: False] := Catch[Module[{f = f0, gens = T["gens"], q = T["q"], nc, rem, found = {}, den = 1,
    monos, p, cl, branch, eta, delta, spec, vpf, places, g, rho, ell, tr, tau, em, monosHere, cands, coeffs, rhoPart, subT, tt, ev, xi, lamS, val, K, v,
    kern, lp, rp, piece, r0, cP, dl, tauNc},
  nc = NCo[T];
  rem = TZero[T];
  Do[den = PolynomialLCM[den, Denominator[Can[f[[i]]]]], {i, nc}];
  If[FreeQ[den, Alternatives @@ gens], Throw[{f, rem, found}, "rs"]];
  monos = If[special && sources =!= {}, Monomials[T, sources], {}];
  Do[p = fac[[1]];
    If[FreeQ[p, Alternatives @@ gens], Continue[]];
    cl = QuietCheck[ClassifyPrime[T, p], $Failed];
    If[cl === $Failed, Continue[]];
    {branch, eta, delta, spec} = cl;
    If[spec || branch, Continue[]];
    vpf = VP[f, p, T, branch];
    If[vpf > -delta, Continue[]];
    places = LinearPlaces[T, p];
    If[places === None || (vpf < -delta && (q =!= None || Length[places] != 1 || Expand[places[[1, 3]] - p] =!= 0)), Continue[]];
    Do[{g, rho, ell} = pl;
      If[vpf < -delta,
        tr = CanonicalResidueField[T, f, p, g, delta, False];
        If[tr === None, Continue[]];
        tau = TScalar[T, Can[tr]],
        tau = ResidueAt[T, f, g, rho, ell]];
      If[IsConst[T, tau, g, rho], Continue[]];
      em = TauEmax[T, tau];
      monosHere = If[em <= 2, monos, If[special && sources =!= {}, Monomials[T, sources, em], {}]];
      cands = EiCandidates[T, monosHere, g, rho, ell];
      (* the whole residue, else term by term: the terms in the exponential span go
         to Ei kernels, the others (rho) to the remainder *)
      coeffs = If[cands === {}, None, DecompResidue[T, tau, cands]];
      rhoPart = None;
      If[coeffs === None && AllTrue[Rest[tau], Can[#] === 0 &],
        coeffs = ConstantArray[0, Length[cands] + 1];
        rhoPart = 0;
        Do[tt = Prepend[ConstantArray[0, Length[tau] - 1], term];
          subT = If[cands === {}, None, DecompResidue[T, tt, cands]];
          If[subT === None, rhoPart += term, coeffs = coeffs + subT],
          {term, AdditiveTerms[Can[tau[[1]]], gens]}];
        If[AllTrue[Rest[coeffs], ExactZeroQ], coeffs = None]];           (* nothing in the span: the polar part below *)
      If[coeffs =!= None,
        Do[If[ExactZeroQ[coeffs[[j + 1]]], Continue[]];
          {ev, xi, lamS, val, K} = cands[[j, 1 ;; 5]];
          v = Can[lamS - val];
          kern = MakeKernel["gamma", K, <|"s" -> 0, "v" -> v, "omega" -> CanT[TPad[T, TDiv[T, xi, TScalar[T, v]]]], "cw" -> 1,
            "shift" -> -val, "J" -> ev, "place" -> {g, rho}|>];
          AppendTo[found, {coeffs[[j + 1]], kern}];
          f = Padd[f, Pscale[-coeffs[[j + 1]], K]],
          {j, Length[cands]}];
        If[rhoPart === None || Can[rhoPart] === 0, Continue[]];
        (* the part of the residue outside the span: rho~ D ell / ell (residue rho) *)
        lp = TScalar[T, ell];
        rp = Can[rhoPart];
        piece = CanT[TPad[T, TMul[T, TScalar[T, rp], TDiv[T, TowerD[T, lp], lp]]]];
        f = Padd[f, Pscale[-1, piece]];
        rem = Padd[rem, piece];
        If[verbose, Print["  (", p, ") at ", g, " = ", rho, ": residue part ", rp, " carried by the remainder"]];
        Continue[]];
      (* the residue goes to the remainder: at a critical pole, as the polar part
         c_P / ell, c_P = (f ell)|_{g = rho}; at a deeper pole, as tau~ D ell / ell *)
      lp = TScalar[T, ell];
      r0 = Total[Select[If[Head[#] === Plus, List @@ #, {#}] &[Expand[Can[tau[[1]]]]], FreeQ[#, Alternatives @@ gens] &]];   (* the constant part stays *)
      If[vpf == -delta,
        cP = Can[Quiet[Can[# ell] /. g -> rho]] & /@ TPad[T, f];
        piece = Can[#/ell] & /@ cP;
        If[r0 =!= 0,
          dl = TPad[T, TowerD[T, lp]];
          piece = Padd[piece, Pscale[-r0, Can[Quiet[Can[#] /. g -> rho]/ell] & /@ dl]]],
        tauNc = Prepend[Rest[tau], Can[tau[[1]] - r0]];
        piece = TMul[T, TPad[T, tauNc], TDiv[T, TowerD[T, lp], lp]]];
      piece = CanT[TPad[T, piece]];
      f = Padd[f, Pscale[-1, piece]];
      rem = Padd[rem, piece];
      If[verbose, Print["  (", p, ") at ", g, " = ", rho, ": residue ", tau, " carried by the remainder"]],
      {pl, places}],
    {fac, FactorsOf[den]}];
  {CanRT[f], CanRT[rem], found}], "rs"];

(* ReducedColumn[T, col]: (R1) for a remainder column: at most a critical pole at
   every normal prime of its denominator *)
ReducedColumn[T_, col_] := Module[{gens = T["gens"], den = 1, cl},
  Do[den = PolynomialLCM[den, Denominator[Can[z]]], {z, col}];
  If[FreeQ[den, Alternatives @@ gens], Return[True]];
  Do[cl = QuietCheck[ClassifyPrime[T, fac[[1]]], $Failed];
    If[cl === $Failed, Return[False, Module]];
    If[! cl[[4]] && VP[TPad[T, col], fac[[1]], T, cl[[1]]] < -cl[[3]], Return[False, Module]],
    {fac, FactorsOf[den]}];
  True];

(* AdditiveTerms[c, gens]: the additive terms of a rational function: its partial
   fractions in the main generator of its denominator, the polynomial part split
   into terms *)
AdditiveTerms[c0_, gens_] := Module[{c = Can[c0], n, d, main, g, parts, out = {}, pn, pd, terms},
  terms[e_] := If[Head[e] === Plus, List @@ e, {e}];
  {n, d} = {Numerator[c], Denominator[c]};
  main = Select[gens, ! FreeQ[d, #] &];
  If[main === {}, Return[#/d & /@ terms[Expand[n]]]];
  g = Last[main];
  parts = QuietCheck[terms[Apart[c, g]], $Failed];
  If[parts === $Failed, Return[{c}]];
  Do[{pn, pd} = {Numerator[#], Denominator[#]} &[Can[pt]];
    out = Join[out, #/pd & /@ terms[Expand[pn]]], {pt, parts}];
  out];

(* RemainderBasis[T, rem, order]: the remainder columns: first the additive terms
   of rem itself, then, per coordinate, the monomials up to the numerator degrees
   of rem over the denominator of rem with every normal prime reduced to its
   critical exponent (special primes kept), simplest first *)
RemainderBasis[T_, rem_, order_: "terms"] := Module[{gens = T["gens"], nc = NCo[T], cols = {}, termCols = {}, c, col, n, d, DR, cl, spec, delta, degs, monos, mo},
  Do[c = Can[rem[[i]]];
    If[c === 0, Continue[]];
    Do[col = ConstantArray[0, nc]; col[[i]] = Can[term];
      If[ReducedColumn[T, col], AppendTo[cols, col]; AppendTo[termCols, col]],     (* (R1): no deep pole at a normal prime *)
      {term, AdditiveTerms[c, gens]}],
    {i, nc}];
  Do[c = Can[rem[[i]]];
    If[c === 0, Continue[]];
    {n, d} = {Numerator[c], Denominator[c]};
    DR = 1;
    If[! FreeQ[d, Alternatives @@ gens],
      Do[cl = QuietCheck[ClassifyPrime[T, fac[[1]]], $Failed];
        {spec, delta} = If[cl === $Failed, {True, fac[[2]]}, {cl[[4]], cl[[3]]}];
        DR *= fac[[1]]^If[spec, fac[[2]], Min[fac[[2]], delta]],
        {fac, FactorsOf[d]}]];
    degs = Table[If[! FreeQ[n, g], Exponent[n, g], 0] + If[! FreeQ[d, g], Exponent[d, g], 0], {g, gens}];   (* room for the Hermite-reduced numerator *)
    monos = SortBy[Tuples[Range[0, #] & /@ degs], {Total[#], #} &];
    Do[mo = Times @@ (gens^a);
      col = ConstantArray[0, nc]; col[[i]] = Can[mo/DR];
      AppendTo[cols, col], {a, monos}],
    {i, nc}];
  If[order === "monomials",                      (* the reduced monomial basis first, simplest first *)
    cols = Join[Select[cols, ! MemberQ[termCols, #] &], Select[cols, MemberQ[termCols, #] &]]];
  cols];

(* MoveCriticalTerms[f, T, p]: the additive terms of f with a pole of at least the
   critical order at the prime p, moved to the remainder: {f minus them, them} or None *)
MoveCriticalTerms[f_, T_, p_] := Module[{cl, branch, eta, delta, spec, nc = NCo[T], moved, c, col, places},
  cl = QuietCheck[ClassifyPrime[T, p], $Failed];
  If[cl === $Failed, Return[None]];
  {branch, eta, delta, spec} = cl;
  moved = TZero[T];
  Do[c = Can[f[[i]]];
    If[c === 0, Continue[]];
    Do[col = TZero[T]; col[[i]] = Can[term];
      If[! FreeQ[Denominator[Can[term]], Alternatives @@ Append[Select[T["gens"], ! FreeQ[p, #] &], Unique[]]] &&
          VP[col, p, T, branch] <= -delta,
        places = If[branch, None, LinearPlaces[T, p]];
        If[places =!= None && places =!= {} && AllTrue[places, IsConst[T, ResidueAt[T, col, #[[1]], #[[2]], #[[3]]], #[[1]], #[[2]]] &],
          Continue[]];                     (* a constant residue: a logarithm can carry it *)
        moved = Padd[moved, col]],
      {term, AdditiveTerms[c, T["gens"]]}],
    {i, nc}];
  If[AllZeroQ[moved], Return[None]];
  {CanT[Padd[f, Pscale[-1, moved]]], CanT[moved]}];

(* PartialSolve[key, kernels, split, retry, verbose, order]: Steps (P3)-(P4): the
   Part II ansatz with the kernel columns and the remainder columns, row-reduced
   with the ansatz and kernel columns first and the remainder columns last; free
   unknowns 0.  {I, kernel coefficients, remainder} or None *)
PartialSolve[key_, kernels_, split_, retry_, verbose_, order_: "terms"] := Module[{A = $analyses[key], T, f, Y, gens, q, m, nc, detLogs, denv,
    unitsBase, rem, rootLogs, unkLogs, sunits, bounds, exps, proved, monos, css, units, gammas, betas, kappas, Rcols, rhos, unks, extra,
    sub, neq, V, rat, y, surf, lrTerms, L, I0, coeffs, r},
  T = A["T"]; f = A["f"]; Y = A["Y"]; {gens, q, m} = {T["gens"], T["q"], T["m"]};
  nc = NCo[T];
  {detLogs, denv, unitsBase, rem, rootLogs} = Lookup[A, {"detLogs", "denv", "unitsBase", "rem", "rootLogs"}];
  {unkLogs, sunits} = AnalysisSpecials[key, split, False];
  {bounds, exps, proved} = ProxyBounds[key, unkLogs, kernels, retry, False];
  Do[denv *= unkLogs[[i, 1]]^exps[[i]], {i, Length[unkLogs]}];
  monos = Tuples[Range[0, #] & /@ bounds];
  css = Table[cc[i - 1] @@@ monos, {i, nc}];
  units = Join[sunits, unitsBase];
  gammas = Table[gamma[i], {i, Length[units]}];
  betas = Table[beta[i], {i, Length[unkLogs]}];
  kappas = Table[kappa[i], {i, Length[kernels]}];
  Rcols = RemainderBasis[T, rem, order];
  rhos = Table[rhoR[i], {i, Length[Rcols]}];
  unks = Join[Flatten[css], gammas, betas, kappas, rhos];
  extra = Join[Table[{kappas[[j]], TPad[T, kernels[[j]]["col"]]}, {j, Length[kernels]}],
    Table[{rhos[[j]], TPad[T, Rcols[[j]]]}, {j, Length[Rcols]}]];
  {sub, neq} = SolveKernelSystem[T, rem, denv, units, unkLogs, css, monos, gammas, betas, unks, extra];
  If[sub === None, Return[None]];                    (* not even the remainder columns reach f *)
  V = PadRight[Table[Total[MapThread[#1 Times @@ (gens^#2) &, {css[[i]], monos}]]/denv, {i, nc}], T["n"]];
  rat = Can /@ (V /. sub);
  y = If[q =!= None, q^(1/m), None];
  surf[u_] := If[q === None, u[[1]], Sum[TPad[T, u][[i + 1]] y^i/T["E"][[i + 1]], {i, 0, T["n"] - 1}]];
  lrTerms = Join[{#[[1]], ToY[T, #[[2]], Y]} & /@ detLogs,
    Table[{gammas[[i]] /. sub, ToY[T, units[[i, 1]], Y]}, {i, Length[units]}],
    Table[{betas[[i]] /. sub, unkLogs[[i, 1]]}, {i, Length[unkLogs]}]];
  lrTerms = Select[lrTerms, ! ExactZeroQ[#[[1]]] &];
  {L, rat} = If[lrTerms =!= {}, LogToReal[lrTerms, rat, f, T, Y, verbose], {0, rat}];
  I0 = surf[rat] + If[q === None, L, L /. Y -> y] + Total[rootLogs];
  coeffs = kappas /. sub;
  r = TZero[T];
  Do[If[! ExactZeroQ[rhos[[j]] /. sub], r = Padd[r, Pscale[rhos[[j]] /. sub, TPad[T, Rcols[[j]]]]]], {j, Length[Rcols]}];
  If[verbose, Print["  partial: bounds ", bounds, ", ", Length[unks], " unknowns (", Length[Rcols], " remainder columns), ",
    Count[rhos /. sub, _?(! ExactZeroQ[#] &)], " used"]];
  {I0, coeffs, CanT[r]}];

(* RemainderSize[T, r]: {number of additive terms, total degree, operations} of a
   remainder: the order in which partial answers are compared *)
RemainderSize[T_, r_] := Module[{terms = 0, deg = 0, ops = 0, z, n, d, syms, td},
  td[p_, vs_] := If[p === 0, 0, Max[Total /@ Keys[CoefficientRules[p, vs]]]];
  Do[z = Can[z0];
    If[z === 0, Continue[]];
    {n, d} = {Numerator[z], Denominator[z]};
    terms += Length[If[Head[#] === Plus, List @@ #, {#}] &[Expand[n]]];
    syms = Select[T["gens"], ! FreeQ[z, #] &];
    deg += If[syms =!= {}, td[n, syms] + td[d, syms], 0];
    ops += CountOps[z],
    {z0, r}];
  {terms, deg, ops}];

PartialFinish[res0_, T_, fCert_, special_] := Module[{res = res0},
  res["remainder"] = CanT[res["remainder"]];
  res["remainderCertified"] = fCert;
  res];

Options[PartialIntegrate] = {"Special" -> True, "Sources" -> None, "Back" -> None, "Y" -> None, "Verbose" -> False};

(* PartialIntegrate[f, T, opts]: Algorithm S10.  A complete answer if Algorithm S7
   finds one; otherwise a partial result f = D(I) + r with a reduced remainder
   (Section 8.3).  "Special" -> False gives the elementary partial answer. *)
PartialIntegrate[f00_List, T_Association, OptionsPattern[]] := WithPart2Memo[Module[{special = TrueQ[OptionValue["Special"]], sources = OptionValue["Sources"],
    back = OptionValue["Back"], Y = OptionValue["Y"], verbose = TrueQ[OptionValue["Verbose"]], f, full, fCert, f1, rem, terms, key, f2, th, moved,
    kernels, best = None, out, I0, coeffs, r, used, S, h, out2, sz, doneQ = False},
  f = CanT[TPad[T, f00]];
  full = If[special, PMStage["= ParallelIntegrateSpecial (bracket)", ParallelIntegrateSpecial[f, T, "Sources" -> sources, "Back" -> back, "Y" -> Y, "Verbose" -> verbose, "Strict" -> True]],
    PMStage["Part2 (P0)", Part2[f, T]]];
  If[ExprQ[full] || SpecialResultQ[full], Return[full]];
  fCert = If[special,
    If[ListQ[full] && full =!= {} && MemberQ[{"not elementary", "not in class"}, full[[1]]], full, None],
    If[Certified[full], full, None]];
  If[sources === None && special, sources = ExpSources[T, back, Y, False]];
  If[sources === None, sources = {}];
  (* (P1) residues; if an Ei kernel of the exponential span leaves a new
     non-constant residue elsewhere, the residues are carried by the remainder instead *)
  {f1, rem, terms} = PMStage["ResidueSplits", ResidueSplits[f, T, sources, special, verbose]];
  key = If[! AllZeroQ[f1], PMStage["SPAnalyse", SPAnalyse[f1, T, verbose]], None];
  If[special && terms =!= {} && ListQ[key] && key[[1]] === "not elementary" && ! StringQ[key[[2]]],
    {f1, rem, terms} = PMStage["ResidueSplits", ResidueSplits[f, T, {}, False, verbose]];
    key = If[! AllZeroQ[f1], PMStage["SPAnalyse", SPAnalyse[f1, T, verbose]], None]];
  (* (P2) third kind for a certified non-torsion divisor *)
  If[ListQ[key] && key[[1]] === "not elementary" && StringQ[key[[2]]] && StringContainsQ[key[[2]], "not torsion"],
    {f2, th} = ThirdKind[f1, T, verbose];
    If[th =!= {},
      If[special, terms = Join[terms, th],
        Do[rem = Padd[rem, Pscale[t[[1]], t[[2]]["col"]]], {t, th}]];
      f1 = f2;
      key = If[! AllZeroQ[f1], SPAnalyse[f1, T, verbose], None]]];
  Do[                                                          (* (P2') a prime where the analysis stops *)
    If[! (ListQ[key] && key[[1]] === "not elementary" && ! StringQ[key[[2]]]), Break[]];
    moved = MoveCriticalTerms[f1, T, key[[2]]];
    If[moved === None, Break[]];
    {f1, h} = moved;
    rem = Padd[rem, h];
    If[verbose, Print["  partial: the terms of f with a critical pole at (", key[[2]], ") are carried by the remainder"]];
    key = If[! AllZeroQ[f1], SPAnalyse[f1, T, verbose], None],
    {4}];
  If[key === None, Return[PartialFinish[PartialResult[0, terms, rem], T, fCert, special]]];
  If[ListQ[key],
    If[verbose, Print["  partial: the remainder analysis stops (", StringTake[ToStr[key], UpTo[80]], "); r = f"]];
    Return[PartialFinish[PartialResult[0, {}, f], T, fCert, special]]];
  kernels = If[special, Join[PMStage["GammaCandidates", GammaCandidates[f1, T, sources, False]], PMStage["EllipticColumns", EllipticColumns[T, False]]], {}];
  Do[Do[                                                       (* both splittings and both column orders; the simplest remainder wins *)
      out = PMStage["PartialSolve", PartialSolve[key, kernels, split, 0, verbose, order]];
      If[out === None, Continue[]];
      {I0, coeffs, r} = out;
      used = Select[Transpose[{coeffs, kernels}], ! ExactZeroQ[#[[1]]] &];
      If[used =!= {},                                          (* strict: special terms only if necessary *)
        S = SpecialPart[SpecialResult[0, used], T];
        h = CanT[Padd[f1, Pscale[-1, r]]];                      (* the part that I integrates *)
        If[PMStage["CertifyNonelementary", CertifyNonelementary[S, T]] === None && PMStage["CertifyNonelementary", CertifyNonelementary[h, T]] === None,
          out2 = PMStage["PartialSolve", PartialSolve[key, {}, split, 0, verbose, order]];
          If[out2 =!= None, {I0, coeffs, r} = out2; used = {}]]];
      sz = RemainderSize[T, r];
      If[best === None || LexLess[sz, best[[4]]], best = {I0, used, r, sz}];
      If[sz[[1]] == 0, Break[]],
      {order, {"terms", "monomials"}}];
    If[best =!= None && best[[4, 1]] == 0, Break[]],
    {split, {False, True}}];
  If[best === None, Return[PartialFinish[PartialResult[0, {}, f], T, fCert, special]]];
  {I0, used, r} = best[[1 ;; 3]];
  PartialFinish[PartialResult[I0, Join[terms, used], Padd[rem, r]], T, fCert, special]]];

Options[IntegrateSurfacePartial] = {"Special" -> True, "Verbose" -> False, "Tower" -> None,
  "Samples" -> None, "Additive" -> True};

(* IntegrateSurfacePartial[f, x, opts]: IntegrateSurfaceSpecial with a partial
   answer I + Inactive[Integrate][r, x] when no complete one is found;
   {answer, verified} or a status *)
IntegrateSurfacePartial[integrand_, x_Symbol, OptionsPattern[]] := TimeConstrained[Module[{special = TrueQ[OptionValue["Special"]], verbose = TrueQ[OptionValue["Verbose"]],
    tower = OptionValue["Tower"], samples = OptionValue["Samples"], integrandN = integrand /. PMExp -> Exp, bt, T, fp, back, Ysym, res, x0, yExpr, surf, total, anti, fixed, r, ok, d, blocks, sp},
  (* the additive fast path: integrate blocks of terms with pairwise-disjoint
     generator families on their own towers and sum.  Skipped when the caller
     pinned a tower -- honouring a tower for the whole integrand and splitting
     it are contradictory -- and when it declines, the joint path below runs
     unchanged. *)
  If[tower === None && TrueQ[OptionValue["Additive"]],
    blocks = Quiet[AdditiveBlocks[integrand, x]];
    If[ListQ[blocks] && Length[blocks] >= 2,
      sp = SplitIntegratePartial[blocks, integrandN, x, special, verbose, samples];
      If[sp =!= None, Return[sp]]]];
  If[tower === None,
    bt = PMStage["BuildTower", Catch[Quiet[BuildTower[integrand, x, "StructureTheorem" -> True, "Verbose" -> verbose]], "build"]];
    If[! ListQ[bt] || Length[bt] != 4, Return[{"failed", "tower construction failed", ToStr[integrand]}]];
    {T, fp, back} = bt[[1 ;; 3]],
    {T, fp, back} = tower];
  Ysym = BackY[back, T["gens"]];
  res = PMStage["= PartialIntegrate (bracket)", PartialIntegrate[fp, T, "Special" -> special, "Back" -> back, "Y" -> Ysym, "Verbose" -> verbose]];
  If[samples === None,
    x0 = SamplePoint[integrandN, x];
    samples = {If[x0 =!= None, x0, 7/5], 9/4, 13/5}];
  (* a complete answer: Python reruns integrate_surface_special on the same tower,
     which returns the same answer (the pipeline is deterministic, and a special answer
     of PartialIntegrate is certified); it is presented here without the rerun *)
  If[! PartialResultQ[res],
    Return[If[special, PMStage["PresentSpecial", PresentSpecial[res, integrandN, x, T, back, samples, $StrictSP, False, False]],
      If[ExprQ[res], {SubstituteBack[res, back], True}, res]]]];
  yExpr = If[T["q"] =!= None, T["q"]^(1/T["m"]), None];
  surf[u0_] := With[{u = TPad[T, u0]},
    SubstituteBack[u[[1]] + If[T["q"] =!= None, Sum[u[[i + 1]] yExpr^i/T["E"][[i + 1]], {i, 1, T["n"] - 1}], 0], back]];
  total = SubstituteBack[res["elementary"], back];
  Do[anti = SubstituteBack[KernelAntiderivative[t[[2]]], back];
    fixed = FixBranch[t[[1]] anti, t[[1]] surf[t[[2]]["col"]], x, samples, t[[2]]];
    total += If[fixed =!= None, fixed, t[[1]] anti],
    {t, res["terms"]}];
  total = PMStage["Present (partial)", HalfAngleFold[Present[total], back]];
  r = PMStage["Simplify (remainder)", Simplify[surf[res["remainder"]]]];
  ok = False;
  d = PMStage["VerifyAnswer (partial)", QuietCheck[D[total, x] + r, $Failed]];
  If[d =!= $Failed,
    ok = AllTrue[samples, With[{a = Quiet[Num[d, x, #]], b = Quiet[Num[integrandN, x, #]]},
      NumericQ[a] && NumericQ[b] && Abs[a - b] < 10^-15 (1 + Abs[b])] &]];
  {total + Inactive[Integrate][r, x], ok}],
  $SpecialTimeBudget, {"failed", "time budget exceeded"}];


End[];

(* The Integrate`ParallelMixedSpecial method symbol is a C builtin (see
   src/calculus/integrate.c: builtin_integrate_pms) that lazily loads this file
   on first use and delegates to ParallelMixed`Private`IntegrateSurfaceSpecial,
   so the package load is paid only by sessions that actually reach the method. *)
