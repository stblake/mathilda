(* ::Package:: *)

(* ========================================================================== *)
(* Common.m -- the failure protocol and small shared predicates.              *)
(*                                                                            *)
(* A stage that cannot continue THROWS a tragerFail[class, reason, info]       *)
(* object with tag "TragerFail"; the public entry point catches it and turns   *)
(* it into the list it returns.  There are five classes:                      *)
(*                                                                            *)
(*   "NonElementary" -- a proof that the integral is not elementary (the      *)
(*                      reason names the certificate);                        *)
(*   "Incomplete"    -- the integrand is outside the class the procedure      *)
(*                      decides (nested radicals, a radical compositum, a     *)
(*                      general algebraic function): "implementation          *)
(*                      incomplete", reported to the user;                    *)
(*   "Unimplemented" -- inside that class, but a branch of the algorithm is    *)
(*                      not written yet.  Reported like "Incomplete"; each    *)
(*                      site carries an id that tools/check_bronstein_        *)
(*                      completeness.py ratchets to zero;                     *)
(*   "Decline"       -- a resource bound or a budget gave out: no verdict,    *)
(*                      silent;                                               *)
(*   "Internal"      -- an invariant failed.  This is a bug; the entry point   *)
(*                      still only declines, so a bug can never become a       *)
(*                      wrong answer or a false certificate.                  *)
(* ========================================================================== *)

tragerFail::usage = "tragerFail[class, reason, info] is the inert object thrown (tag \"TragerFail\") when a stage stops: class is \"NonElementary\", \"Decline\" or \"Internal\", reason a short string naming the cause, info an Association of diagnostic data.";

throwNonElementary::usage = "throwNonElementary[reason, info] throws tragerFail[\"NonElementary\", reason, info]: a certificate that the integral is not elementary. Use it only where a theorem of [T] proves non-elementarity.";
throwNonElementary[reason_String, info_Association : <||>] :=
  Throw[tragerFail["NonElementary", reason, info], "TragerFail"];

throwDecline::usage = "throwDecline[reason, info] throws tragerFail[\"Decline\", reason, info]: the integral is out of scope or a resource bound was reached; no verdict is claimed.";
throwDecline[reason_String, info_Association : <||>] :=
  Throw[tragerFail["Decline", reason, info], "TragerFail"];

throwIncomplete::usage = "throwIncomplete[reason, info] throws tragerFail[\"Incomplete\", reason, info]: the integrand is outside the class the procedure decides (one simple radical over an elementary tower). Integrate reports it as \"implementation incomplete\".";
throwIncomplete[reason_String, info_Association : <||>] :=
  Throw[tragerFail["Incomplete", reason, info], "TragerFail"];

throwUnimplemented::usage = "throwUnimplemented[site, reason, info] throws tragerFail[\"Unimplemented\", reason, <|\"site\" -> site, ...|>]: the integrand is inside the decided class but this branch of the algorithm is not implemented yet. Integrate reports it as \"implementation incomplete\". site is a stable id (\"File:tag\") counted by the completeness gate, which ratchets the sites to zero.";
throwUnimplemented[site_String, reason_String, info_Association : <||>] :=
  Throw[tragerFail["Unimplemented", reason, Append[info, "site" -> site]], "TragerFail"];

throwInternal::usage = "throwInternal[invariant, info] throws tragerFail[\"Internal\", invariant, info]: an invariant that the mathematics guarantees has failed. The entry point declines on it, so a bug cannot become a wrong answer.";
throwInternal[invariant_String, info_Association : <||>] :=
  Throw[tragerFail["Internal", invariant, info], "TragerFail"];

catchTrager::usage = "catchTrager[expr] evaluates expr and returns its value, or the tragerFail[...] object thrown by any stage inside it.";
SetAttributes[catchTrager, HoldFirst];
catchTrager[expr_] := Catch[expr, "TragerFail"];

tragerFailQ::usage = "tragerFailQ[e] is True when e is a tragerFail[...] object.";
tragerFailQ[e_] := MatchQ[e, _tragerFail];

tragerRetryableQ::usage = "tragerRetryableQ[e] is True when e is a tragerFail[...] object that proves nothing (class Decline, Unimplemented or Incomplete), so that a caller may try another construction; False for a certificate (NonElementary), an Internal error, or a non-failure.";
tragerRetryableQ[tragerFail["Decline" | "Unimplemented" | "Incomplete", _, _]] := True;
tragerRetryableQ[_] := False;

failToList::usage = "failToList[tragerFail[class, reason, info]] is the list returned to Integrate: {\"not elementary\", reason} for a certificate, {\"implementation incomplete\", reason} for an Incomplete or Unimplemented stop, {\"failed\", reason} otherwise (an internal error is reported as {\"failed\", \"internal: \" <> invariant}).";
failToList[tragerFail["NonElementary", reason_, _]] := {"not elementary", reason};
failToList[tragerFail["Incomplete", reason_, _]] := {"implementation incomplete", reason};
failToList[tragerFail["Unimplemented", reason_, _]] := {"implementation incomplete", reason};
failToList[tragerFail["Decline", reason_, _]] := {"failed", reason};
failToList[tragerFail["Internal", reason_, _]] := {"failed", "internal: " <> reason};

(* -------------------------------------------------------------------------- *)
(* Predicates on the base field Q.                                            *)
(* -------------------------------------------------------------------------- *)

rationalNumberQ::usage = "rationalNumberQ[c] is True for an explicit integer or rational number.";
rationalNumberQ[c_] := IntegerQ[c] || Head[c] === Rational;

rationalPolyQ::usage = "rationalPolyQ[p, x] is True when p is a polynomial in x with coefficients in the constant field (after expansion): rationals, and Gaussian rationals under tragerConstantI.";
rationalPolyQ[p_, x_Symbol] :=
  PolynomialQ[p, x] && AllTrue[CoefficientList[Expand[p], x], constantFieldNumberQ];

constantFieldNumberQ::usage = "constantFieldNumberQ[c] is True when c is in the constant field in use: a rational number; a Gaussian rational while tragerConstantI is True (Q(i), Phase 4.1); an explicit algebraic number built from rationals and radicals of rationals while tragerAlgConstants is non-empty (Phase 4.3).";
constantFieldNumberQ[c_] := rationalNumberQ[c] ||
  (TrueQ[tragerConstantI] && MatchQ[c, Complex[(_Integer | _Rational), (_Integer | _Rational)]]) ||
  (tragerAlgConstants =!= {} && NumericQ[c] &&
    FreeQ[c, Pi | E | EulerGamma | Catalan | _Log | _Real | _Root | _AlgebraicNumber] &&
    FreeQ[c, Power[b_, e_] /; ! (rationalNumberQ[b] && rationalNumberQ[e]) && ! (rationalNumberQ[e] && IntegerQ[e])] &&
    (TrueQ[tragerConstantI] || FreeQ[c, _Complex]));

rationalFunctionQ::usage = "rationalFunctionQ[e, x] is True when e is a rational function of x with coefficients in the constant field (constantFieldNumberQ).";
rationalFunctionQ[e_, x_Symbol] := Module[{t = Together[e]},
  rationalPolyQ[Numerator[t], x] && rationalPolyQ[Denominator[t], x]];

algebraicNumberQ::usage = "algebraicNumberQ[c] is True when c is an explicit algebraic number: rationals combined with radicals of rationals by + - * / and rational powers (no symbol, no transcendental constant).";
algebraicNumberQ[c_] := NumericQ[c] && FreeQ[c, Pi | E | EulerGamma | Catalan | _Log | _Exp | _Real | _Complex] &&
  FreeQ[c, Power[b_, e_] /; !rationalNumberQ[e]];

algebraicRationalFunctionQ::usage = "algebraicRationalFunctionQ[e, x] is True when e is a rational function of x whose coefficients are algebraic numbers (algebraicNumberQ): a rational function over a number field Q(theta).";
algebraicRationalFunctionQ[e_, x_Symbol] := Module[{t = Together[e], num, den},
  num = Expand[Numerator[t]]; den = Expand[Denominator[t]];
  PolynomialQ[num, x] && PolynomialQ[den, x] &&
    AllTrue[Join[CoefficientList[num, x], CoefficientList[den, x]], algebraicNumberQ]];

rationalInXYQ::usage = "rationalInXYQ[e, x, y] is True when e is a rational function of the two symbols x, y with coefficients in the constant field (constantFieldNumberQ).";
rationalInXYQ[e_, x_Symbol, y_Symbol] := Module[{t = Together[e], num, den},
  num = Expand[Numerator[t]]; den = Expand[Denominator[t]];
  PolynomialQ[num, {x, y}] && PolynomialQ[den, {x, y}] &&
    AllTrue[Flatten[{CoefficientList[num, {x, y}], CoefficientList[den, {x, y}]}],
      constantFieldNumberQ]];

rdeg::usage = "rdeg[r, x] is the degree of the rational function r in x (degree of numerator minus degree of denominator); -Infinity for 0.";
rdeg[r_, x_Symbol] := Module[{t = Together[r]},
  If[t === 0, -Infinity,
    Exponent[Numerator[t], x] - Exponent[Denominator[t], x]]];

tragerOption::usage = "tragerOption[name] is the default value of an internal setting. None of them limits the problems the method decides: searches (shift points, good primes, torsion orders, residue fields) run as long as the problem needs. \"KauersMaxPower\" (largest ideal power tried by Kauers' heuristic, a heuristic method by definition), \"RealFormMaxLeaves\" and \"RealFormTime\" (size and time limits of the real-form presentation; the verified answer is returned either way), \"DirectPowerMaxOrder\" (largest torsion order whose generator is found from the ideal power directly rather than by reduced class arithmetic: a choice of method), \"TorsionScreenPrimes\" (extra good primes at which N delta-bar must vanish before the characteristic-0 principality test of N delta is attempted; a nonzero class there proves N delta is not principal), \"NativeKernels\" (True: Hermite forms, Popov forms, determinants/adjugates and Q(theta)[x] division and gcds go through the FLINT` polynomial-matrix kernels first, the package routines remaining the fallback; False: package routines only, the reference for tests and A/B timing).";
tragerOption["KauersMaxPower"] = 16;
tragerOption["RealFormMaxLeaves"] = 1500;
tragerOption["RealFormTime"] = 5;
tragerOption["DirectPowerMaxOrder"] = 12;
tragerOption["NativeKernels"] = True;
tragerOption["TorsionScreenPrimes"] = 3;

tragerConstantI::usage = "tragerConstantI is True while an integrand with Gaussian-rational constants is integrated (constant field Q(i), Bronstein Phase 4.1). Factorisations are then over Q(i) (tragerFactorList), every number field built contains i or is relative over Q(i), and reduction mod p refuses the data (goodPrimes).";
tragerConstantI = False;

tragerFactorList::usage = "tragerFactorList[p] is FactorList[p] over the constant field: over Q(i) (Extension -> I) when tragerConstantI is True, so an irreducible factor is irreducible over the constants actually in use (a factor irreducible over Q but split over Q(i), z^2 + 1, would make Q(i)[t]/(mu) a ring with zero divisors and its zero tests wrong). Over Q otherwise.";
tragerFactorList[p_] := With[{ext = tragerExtension[]},
  If[ext === {}, FactorList[p], FactorList[p, Extension -> ext]]];

tragerTogether::usage = "tragerTogether[e] is Together[e] over the constant field in use: with radical constants (tragerAlgConstants non-empty) Together[e, Extension -> Automatic], since a cancellation that needs Sqrt[2]^2 = 2 is otherwise missed and a combined zero reads as a non-zero sum (Phase 4.3b). The same as rc[fieldQ[], e].";
tragerTogether[e_] := If[tragerAlgConstants === {}, Together[e], Together[e, Extension -> Automatic]];
tragerGCD::usage = "tragerGCD[a, b] is PolynomialGCD[a, b] over the constant field in use (Extension -> Automatic with radical constants, Phase 4.3b).";
tragerGCD[a_, b_] := If[tragerAlgConstants === {}, PolynomialGCD[a, b], PolynomialGCD[a, b, Extension -> Automatic]];

tragerSquareFreeList::usage = "tragerSquareFreeList[p] is FactorSquareFreeList[p] over the constant field in use (Extension -> tragerExtension[] when it is not Q). Without the extension a polynomial with algebraic coefficients is not decomposed at all (FactorSquareFreeList[z^2 - Sqrt[2] z^3] is the polynomial itself, multiplicity 1), and Hermite reduction then saw no repeated factor: a false certificate for (x - 2^(1/3))^(-2/3) (Phase 4.3b, G-nf-8).";
tragerSquareFreeList[p_] := With[{ext = tragerExtension[]},
  If[ext === {}, FactorSquareFreeList[p], FactorSquareFreeList[p, Extension -> ext]]];

tragerAlgConstants::usage = "tragerAlgConstants is the list of the algebraic constants of the integrand other than I -- radicals of rational numbers such as Sqrt[2], 2^(1/3) -- while it is integrated (constant field Q(i?, constants), Bronstein Phase 4.3): factorisations are over the field they generate (tragerFactorList), residue fields contain them, and reduction mod p refuses them (gfMod).";
tragerAlgConstants = {};

tragerExtension::usage = "tragerExtension[] is the list of algebraic generators of the constant field in use, for FactorList's Extension option: I when tragerConstantI is True, then tragerAlgConstants.";
tragerExtension[] := Join[If[TrueQ[tragerConstantI], {I}, {}], tragerAlgConstants];

tragerConstantsOf::usage = "tragerConstantsOf[f] is the list of the radical constants of f: every Power[r, q] with r a rational number and q a non-integer rational (Sqrt[2], 2^(1/3), (-3)^(1/4)), DeleteDuplicates'd.";
tragerConstantsOf[f_] := DeleteDuplicates[Cases[{f}, p : Power[r_, q_Rational] /; rationalNumberQ[r] :> p, Infinity]];

nfDemoteCertificate::usage = "nfDemoteCertificate[r] is the failure r unchanged, except that a NonElementary certificate becomes Unimplemented (site NF:certificate) while radical constants are in use (tragerAlgConstants non-empty, Phase 4.3): the base-field arithmetic over Q is exact for Gaussian rationals but not for radicals (1/(1 + Sqrt[2]) is not rationalised), so a zero can read as nonzero and a verdict is not yet a proof. Answers are unaffected: they are verified.";
nfDemoteCertificate[r_] := If[tragerAlgConstants =!= {} && tragerFailQ[r] && r[[1]] === "NonElementary",
  tragerFail["Unimplemented", "non-elementarity over a number-field constant field is not yet certified (" <> r[[2]] <> ")",
    <|"site" -> "NF:certificate"|>], r];

