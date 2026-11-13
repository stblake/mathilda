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

rationalPolyQ::usage = "rationalPolyQ[p, x] is True when p is a polynomial in x with rational coefficients (after expansion).";
rationalPolyQ[p_, x_Symbol] :=
  PolynomialQ[p, x] && AllTrue[CoefficientList[Expand[p], x], rationalNumberQ];

rationalFunctionQ::usage = "rationalFunctionQ[e, x] is True when e is a rational function of x with rational coefficients.";
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

rationalInXYQ::usage = "rationalInXYQ[e, x, y] is True when e is a rational function of the two symbols x, y with rational coefficients.";
rationalInXYQ[e_, x_Symbol, y_Symbol] := Module[{t = Together[e], num, den},
  num = Expand[Numerator[t]]; den = Expand[Denominator[t]];
  PolynomialQ[num, {x, y}] && PolynomialQ[den, {x, y}] &&
    AllTrue[Flatten[{CoefficientList[num, {x, y}], CoefficientList[den, {x, y}]}],
      rationalNumberQ]];

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
