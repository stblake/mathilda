(* ::Package:: *)

(* ========================================================================== *)
(* Common.m -- the failure protocol and small shared predicates.              *)
(*                                                                            *)
(* A stage that cannot continue THROWS a tragerFail[class, reason, info]       *)
(* object with tag "TragerFail"; the public entry point catches it and turns   *)
(* it into the list it returns.  There are three classes:                     *)
(*                                                                            *)
(*   "NonElementary" -- a proof that the integral is not elementary (the      *)
(*                      reason names the certificate);                        *)
(*   "Decline"       -- out of scope, or a bound or heuristic gave out: no     *)
(*                      verdict;                                              *)
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

throwInternal::usage = "throwInternal[invariant, info] throws tragerFail[\"Internal\", invariant, info]: an invariant that the mathematics guarantees has failed. The entry point declines on it, so a bug cannot become a wrong answer.";
throwInternal[invariant_String, info_Association : <||>] :=
  Throw[tragerFail["Internal", invariant, info], "TragerFail"];

catchTrager::usage = "catchTrager[expr] evaluates expr and returns its value, or the tragerFail[...] object thrown by any stage inside it.";
SetAttributes[catchTrager, HoldFirst];
catchTrager[expr_] := Catch[expr, "TragerFail"];

tragerFailQ::usage = "tragerFailQ[e] is True when e is a tragerFail[...] object.";
tragerFailQ[e_] := MatchQ[e, _tragerFail];

failToList::usage = "failToList[tragerFail[class, reason, info]] is the list returned to Integrate: {\"not elementary\", reason} for a certificate, {\"failed\", reason} otherwise (an internal error is reported as {\"failed\", \"internal: \" <> invariant}).";
failToList[tragerFail["NonElementary", reason_, _]] := {"not elementary", reason};
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
