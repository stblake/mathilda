(* ::Package:: *)
(* PARANOIA -- W. Kahan's floating-point test, ported to Mathilda.

   Source: the C translation by Thos Sumner (UCSF, 1985) and David Gay
   (BTL, 1986).  Phase 9 of the 2026-11 numeric-tower review
   (docs/design/numeric_tower_review.md).

   ParanoiaRun[mode] runs the test under a number "mode" and returns
       <|"Mode", "Compare", "Radix", "Precision", "Failure", "Serious",
         "Defect", "Flaw", "Messages" -> {{class, milestone, text}, ...},
         "Notes" -> {{milestone, text}, ...}, "Policy" -> {...}|>
   and Prints a readable report (option "Verbose" -> False suppresses it).

   Modes
     "Machine"        mk[k] = N[k]          (machine doubles)
     {"MPFR", d}      mk[k] = N[k, d]       (fixed-precision MPFR)
     "Mixed"          mk alternates N[k] / N[k, 40] per constant, so the
                      program exercises machine/MPFR contagion.

   Option "Compare" -> "Exact" (default) | "Language".
     PARANOIA's comparisons are C's exact comparisons.  The language's own
     Equal/Less carry a relative tolerance (the last 7 bits of a machine
     number), so under "Exact" every comparison is made on the exact
     values (sign, binary exponent, exact mantissa); a non-real operand
     (Indeterminate, Underflow[], Overflow[], ComplexInfinity) is unordered,
     like a C NaN.  "Language" uses TrueQ[a == b], TrueQ[a < b], ... and
     shows what a program written naively in the language would see.

   Fidelity rules of the port
     * Every FLOAT constant is made by mk[k]; every other value arises by
       arithmetic on them.  C's left-to-right association is reproduced by
       explicit parentheses, so no multi-argument Plus/Times is ever formed
       whose evaluation order the evaluator might choose.
     * FLOOR returns a FLOAT: fl[x] = x - (x - Floor[x]) keeps the mode's
       number type (Floor[x] alone is an Integer).
     * Milestones 0-100 are a line-by-line port.  Milestones 110-220
       (underflow/overflow) are ADAPTED: Mathilda has no subnormals, a
       machine result outside [DBL_MIN, DBL_MAX] continues as a 53-bit MPFR
       number, and Underflow[]/Overflow[] appear only past MPFR's exponent
       range (~2^(+-2^30)).  The thresholds are found by exponent bisection
       on Radix^k and PARANOIA's own checks are then applied to them, for
       two tiers: the machine tier (machine type boundary) and the absolute
       tier (MPFR exponent range).  Outcomes that are documented policy are
       recorded under "Policy", not as FAILURE/SERIOUS. *)

BeginPackage["Paranoia`"];

ParanoiaRun::usage = "ParanoiaRun[mode] runs Kahan's PARANOIA under mode \"Machine\", {\"MPFR\", d} or \"Mixed\" and returns an Association of findings.";
ParanoiaGate::usage = "ParanoiaGate[] runs every mode with exact comparisons (gated) and \"Machine\" with language comparisons (informational); a FAILURE or SERIOUS not on the KNOWN list fails the gate.";

Begin["`Private`"];

(* ------------------------------------------------------------------ *)
(* Error bookkeeping                                                  *)
(* ------------------------------------------------------------------ *)

FAIL = 0; SER = 1; DEF = 2; FLW = 3;
clsName = {"FAILURE", "SERIOUS", "DEFECT", "FLAW"};
rounded = 1; chopped = 2; other = 0;
NoTrials = 20;
loopCap = 20000;

fmt[x_] := ToString[x, InputForm];

bad[k_, t_] := (errCnt[[k + 1]] = errCnt[[k + 1]] + 1;
   AppendTo[msgs, {clsName[[k + 1]], ms, t}]);
(* ErrCnt[k]++ without a message (TstPtUf does this). *)
errInc[k_] := (errCnt[[k + 1]] = errCnt[[k + 1]] + 1);
tst[k_, v_, t_] := If[! TrueQ[v], bad[k, t]];
note[t_] := AppendTo[notes, {ms, t}];
policy[t_] := AppendTo[pol, {ms, t}];
pause[] := (ms = ms + 1);

(* ------------------------------------------------------------------ *)
(* Comparisons                                                        *)
(* ------------------------------------------------------------------ *)

(* Exact key {sign, exponent, |mantissa| as an exact rational} so that a
   comparison never builds the exact value of a number like 2^-(2^30)
   (a gigabit denominator).  x = sign * m * 2^e with m in [1/2, 1). *)
meKey[x_, de_] := Module[{me, m},
   If[Head[x] =!= Real, Return[$bad]];
   me = Quiet[MantissaExponent[x, 2]];
   If[! MatchQ[me, {_Real, _Integer}], Return[$bad]];
   m = SetPrecision[me[[1]], Infinity];
   Which[
     m == 0, If[TrueQ[x == 0], {0, 0, 0}, $bad],
     1/2 <= Abs[m] < 1, {Sign[m], me[[2]] + de, Abs[m]},
     True, $bad]];
exv[x_] := Which[
   Head[x] === Real,
     Module[{k},
       k = meKey[x, 0];
       (* MantissaExponent[x, 2] returns a zero mantissa for an MPFR x in
          the top binade (x = 2^(2^30-2) at 30 digits): retry on the
          exactly scaled x/2 and 2x, and validate every answer. *)
       If[k === $bad, k = meKey[x/2, 1]];
       If[k === $bad, k = meKey[2*x, -1]];
       If[k === $bad, $unord, k]],
   Head[x] === Integer || Head[x] === Rational,
     If[x == 0, {0, 0, 0},
       Module[{e = Floor[Log[2., Abs[x]]] + 1, m},
         m = Abs[x]/2^e;
         While[m >= 1, e = e + 1; m = m/2];
         While[m < 1/2, e = e - 1; m = 2 m];
         {Sign[x], e, m}]],
   True, $unord];
cmpv[a_, b_] := Module[{x, y},
   x = exv[a]; y = exv[b];
   Which[x === $unord || y === $unord, $unord,
     x[[1]] =!= y[[1]], Sign[x[[1]] - y[[1]]],
     x[[1]] === 0, 0,
     x[[2]] =!= y[[2]], x[[1]]*Sign[x[[2]] - y[[2]]],
     True, x[[1]]*Sign[x[[3]] - y[[3]]]]];

eq[a_, b_] := If[$cmp === "Language", TrueQ[a == b], cmpv[a, b] === 0];
ne[a_, b_] := ! eq[a, b];
lt[a_, b_] := If[$cmp === "Language", TrueQ[a < b], cmpv[a, b] === -1];
gt[a_, b_] := If[$cmp === "Language", TrueQ[a > b], cmpv[a, b] === 1];
le[a_, b_] := If[$cmp === "Language", TrueQ[a <= b], MemberQ[{-1, 0}, cmpv[a, b]]];
ge[a_, b_] := If[$cmp === "Language", TrueQ[a >= b], MemberQ[{1, 0}, cmpv[a, b]]];

(* ------------------------------------------------------------------ *)
(* Loops (guarded) and the C library                                  *)
(* ------------------------------------------------------------------ *)

SetAttributes[doWhile, HoldAll];
doWhile[body_, cond_] := Module[{n = 0, go = True},
   While[go, body; n = n + 1;
     If[! TrueQ[cond], go = False];
     If[n > loopCap, go = False;
       bad[FAIL, "loop did not terminate within " <> ToString[loopCap] <> " iterations"]]]];

SetAttributes[whileG, HoldAll];
whileG[cond_, body_] := Module[{n = 0, go = True},
   While[go && TrueQ[cond], body; n = n + 1;
     If[n > loopCap, go = False;
       bad[FAIL, "loop did not terminate within " <> ToString[loopCap] <> " iterations"]]]];

fl[x_] := x - (x - Floor[x]);
sgn[x_] := If[ge[x, Zero], One, MinusOne];
pw[x_, y_] := Quiet[Power[x, y]];
toF[i_] := (One + i) - One;            (* (FLOAT) I *)

(* ------------------------------------------------------------------ *)
(* Helpers of the C program                                           *)
(* ------------------------------------------------------------------ *)

rnd[] := Module[{x, y},
   x = Random1 + Random9; y = x*x; y = y*y; x = x*y; y = x - fl[x];
   Random1 = y + (x*Rnd5); Random1];

sqXMinX[errKind_] := Module[{xa, xb},
   xb = X*BInvrse; xa = X - xb;
   SqEr = ((Sqrt[X*X] - xb) - xa)/OneUlp;
   If[ne[SqEr, Zero],
     If[lt[SqEr, MinSqEr], MinSqEr = SqEr];
     If[gt[SqEr, MaxSqEr], MaxSqEr = SqEr];
     J = J + One;
     bad[errKind, "sqrt(X*X) - X = " <> fmt[SqEr] <> " ulps for X = " <> fmt[X]]]];

newD[] := (
   X = Z1*Q;
   X = (fl[Half - (X/Radix)]*Radix) + X;
   Q = ((Q - (X*Z))/Radix) + ((X*X)*(dD/Radix));
   Z = Z - ((Two*X)*dD);
   If[le[Z, Zero], Z = -Z; Z1 = -Z1];
   dD = Radix*dD);

sr3750[] := If[! (lt[X - Radix, Z2 - Radix] || gt[X - Z2, W - Z2]),
   iI = iI + 1;
   X2 = Sqrt[X*dD];
   Y2 = (X2 - Z2) - (Y - Z2);
   X2 = X8/(Y - Half);
   X2 = X2 - ((Half*X2)*X2);
   SqEr = (Y2 + Half) + (Half - X2);
   If[lt[SqEr, MinSqEr], MinSqEr = SqEr];
   SqEr = Y2 - X2;
   If[gt[SqEr, MaxSqEr], MaxSqEr = SqEr]];

isYeqX[] := If[ne[Y, X],
   If[nN <= 0,
     If[eq[Z, Zero] && le[Q, Zero],
       note["WARNING: computing (" <> fmt[Z] <> ")^(" <> fmt[Q] <> ") yielded " <>
         fmt[Y] <> "; which compared unequal to correct " <> fmt[X]],
       bad[DEF, "computing (" <> fmt[Z] <> ")^(" <> fmt[Q] <> ") yielded " <> fmt[Y] <>
         "; which compared unequal to correct " <> fmt[X] <> "; they differ by " <> fmt[Y - X]]]];
   nN = nN + 1];

sr3980[] := Module[{done = False},
   doWhile[
     Q = toF[iI]; Y = pw[Z, Q]; isYeqX[];
     iI = iI + 1;
     If[iI > mM, done = True, X = Z*X],
     ! done && lt[X, W]]];

printIfNPositive[] := If[nN > 0,
   note["Similar discrepancies have occurred " <> ToString[nN] <> " times."]];

tstPtUf[] := (
   nN = 0;
   If[ne[Z, Zero],
     Q9 = (Z + Z)/Z;
     If[lt[Abs[Q9 - Two], Radix*U2],
       note["(Z + Z) / Z = " <> fmt[Q9] <> " for Z = " <> fmt[Z] <> ": O.K."],
       If[lt[Q9, One] || gt[Q9, Two] || exv[Q9] === $unord,
         nN = 1; errInc[SER];
         AppendTo[msgs, {"SERIOUS", ms, "(Z + Z) / Z = " <> fmt[Q9] <> " for Z = " <> fmt[Z] <>
            ": This is a VERY SERIOUS DEFECT!"}],
         nN = 1; errInc[DEF];
         AppendTo[msgs, {"DEFECT", ms, "(Z + Z) / Z = " <> fmt[Q9] <> " for Z = " <> fmt[Z] <>
            ": This is a DEFECT!"}]]];
     V9 = Z*One; Random1 = V9;
     V9 = One*Z; Random2 = V9;
     V9 = Z/One;
     If[eq[Z, Random1] && eq[Z, Random2] && eq[Z, V9],
       If[nN > 0, pause[]],
       nN = 1;
       bad[DEF, "What prints as Z = " <> fmt[Z] <> " compares different from Z * 1 = " <>
          fmt[Random1] <> ", 1 * Z = " <> fmt[Random2] <> " or Z / 1 = " <> fmt[V9]];
       If[ne[Random2, Random1], errInc[DEF];
         bad[DEF, "Multiplication does not commute! 1 * Z = " <> fmt[Random2] <>
            " differs from Z * 1 = " <> fmt[Random1]]];
       pause[]]]);

(* ------------------------------------------------------------------ *)
(* Constants                                                          *)
(* ------------------------------------------------------------------ *)

makeConstants[] := (
   Zero = mk[0]; One = mk[1]; Two = mk[2]; Three = mk[3]; Four = mk[4];
   Five = mk[5]; Eight = mk[8]; Nine = mk[9]; TwentySeven = mk[27];
   ThirtyTwo = mk[32]; TwoForty = mk[240]; MinusOne = mk[-1];
   Half = mk[1/2]; OneAndHalf = mk[3/2];
   (* the C program's few double literals *)
   Rp01 = mk[1/100];            (* 0.01      in Radix = FLOOR(0.01 + U2/U1) *)
   Rnd5 = mk[5/1000000];        (* 0.000005  in Random()                    *)
   Three3 = mk[3];              (* 3.0       in Random9 = SQRT(3.0)         *)
   Ten = mk[10];                (* 10.       in LOG(10.)                    *)
   Thousandth = mk[1/1000]);    (* 0.001                                    *)

(* ------------------------------------------------------------------ *)
(* Milestones 0 - 30                                                  *)
(* ------------------------------------------------------------------ *)

m7[] := (
   ms = 7;
   tst[FAIL, eq[Zero + Zero, Zero], "0+0 != 0"];
   tst[FAIL, eq[One - One, Zero], "1-1 != 0"];
   tst[FAIL, gt[One, Zero], "1 <= 0"];
   tst[FAIL, eq[One + One, Two], "1+1 != 2"];
   Z = -Zero;
   If[ne[Z, 0],
     errInc[FAIL];
     AppendTo[msgs, {"FAILURE", ms, "Comparison alleges that -0.0 is Non-zero!"}];
     U2 = Thousandth; Radix = 1; tstPtUf[]];
   tst[FAIL, eq[Three, Two + One], "3 != 2+1"];
   tst[FAIL, eq[Four, Three + One], "4 != 3+1"];
   tst[FAIL, eq[Four + (Two*(-Two)), Zero], "4+2*(-2) != 0"];
   tst[FAIL, eq[(Four - Three) - One, Zero], "4-3-1 != 0"];
   tst[FAIL, eq[MinusOne, 0 - One] && eq[MinusOne + One, Zero] &&
      eq[One + MinusOne, Zero] && eq[MinusOne + Abs[One], Zero] &&
      eq[MinusOne + (MinusOne*MinusOne), Zero],
     "-1+1 != 0, (-1)+abs(1) != 0, or -1+(-1)*(-1) != 0"];
   tst[FAIL, eq[(Half + MinusOne) + Half, Zero], "1/2 + (-1) + 1/2 != 0"]);

m10[] := (
   ms = 10;
   tst[FAIL, eq[Nine, Three*Three] && eq[TwentySeven, Nine*Three] &&
      eq[Eight, Four + Four] && eq[ThirtyTwo, Eight*Four] &&
      eq[(((ThirtyTwo - TwentySeven) - Four) - One), Zero],
     "9 != 3*3, 27 != 9*3, 32 != 8*4, or 32-27-4-1 != 0"];
   tst[FAIL, eq[Five, Four + One] && eq[TwoForty, ((Four*Five)*Three)*Four] &&
      eq[(TwoForty/Three) - ((Four*Four)*Five), Zero] &&
      eq[(TwoForty/Four) - ((Five*Three)*Four), Zero] &&
      eq[(TwoForty/Five) - ((Four*Three)*Four), Zero],
     "5 != 4+1, 240/3 != 80, 240/4 != 60, or 240/5 != 48"];
   (* Searching for Radix and Precision. *)
   W = One;
   doWhile[W = W + W; Y = W + One; Z = Y - W; Y = Z - One,
     lt[MinusOne + Abs[Y], Zero]];
   prec = Zero; Y = One;
   doWhile[Radix = W + Y; Y = Y + Y; Radix = Radix - W, eq[Radix, Zero]];
   If[lt[Radix, Two], Radix = One];
   note["Radix = " <> fmt[Radix]];
   If[ne[Radix, 1],
     W = One;
     doWhile[prec = prec + One; W = W*Radix; Y = W + One, eq[Y - W, One]]];
   U1 = One/W; U2 = Radix*U1;
   note["Closest relative separation found is U1 = " <> fmt[U1]];
   E0 = Radix; E1 = U1; E9 = U2; E3 = prec;
   X = Four/Three; Third = X - One; F6 = Half - Third; X = F6 + F6;
   X = Abs[X - Third];
   If[lt[X, U2], X = U2];
   doWhile[U2 = X; Y = (Half*U2) + ((ThirtyTwo*U2)*U2); Y = One + Y; X = Y - One,
     ! (le[U2, X] || le[X, Zero])];
   X = Two/Three; F6 = X - Half; Third = F6 + F6; X = Third - Half;
   X = Abs[X + F6];
   If[lt[X, U1], X = U1];
   doWhile[U1 = X; Y = (Half*U1) + ((ThirtyTwo*U1)*U1); Y = Half - Y; X = Half + Y;
     Y = Half - X; X = Half + Y,
     ! (le[U1, X] || le[X, Zero])];
   If[eq[U1, E1], note["confirms closest relative separation U1."],
     note["gets better closest relative separation U1 = " <> fmt[U1]]];
   W = One/U1; F9 = (Half - U1) + Half;
   Radix = fl[Rp01 + (U2/U1)];
   If[eq[Radix, E0], note["Radix confirmed."],
     note["MYSTERY: recalculated Radix = " <> fmt[Radix]]];
   tst[DEF, le[Radix, Eight + Eight], "Radix is too big: roundoff problems"];
   tst[FLW, eq[Radix, Two] || eq[Radix, 10] || eq[Radix, One],
     "Radix is not as good as 2 or 10"]);

m20[] := (
   ms = 20;
   tst[FAIL, lt[F9 - Half, Half], "(1-U1)-1/2 < 1/2 is FALSE, prog. fails?"];
   X = F9; iI = 1; Y = X - Half; Z = Y - Half;
   tst[FAIL, ne[X, One] || eq[Z, Zero], "Comparison is fuzzy,X=1 but X-1/2-1/2 != 0"];
   X = One + U2; iI = 0;
   ms = 25;
   BMinusU2 = Radix - One; BMinusU2 = (BMinusU2 - U2) + One;
   If[ne[Radix, One],
     X = ((-TwoForty)*Log[U1])/Log[Radix];
     Y = fl[Half + X];
     If[lt[Abs[X - Y]*Four, One], X = Y];
     prec = X/TwoForty;
     Y = fl[Half + prec];
     If[lt[Abs[prec - Y]*TwoForty, Half], prec = Y]];
   If[ne[prec, fl[prec]] || eq[Radix, One],
     note["Precision cannot be characterized by an Integer number of significant digits."]];
   If[eq[Radix, One],
     note["logarithmic encoding has precision characterized solely by U1."],
     note["The number of significant digits of the Radix is " <> fmt[prec]]];
   tst[SER, lt[((U2*Nine)*Nine)*TwoForty, One], "Precision worse than 5 decimal figures"]);

m30[] := (
   ms = 30;
   X = Abs[((((Four/Three) - One) - (One/Four))*Three) - (One/Four)];
   doWhile[Z2 = X; X = (One + ((Half*Z2) + ((ThirtyTwo*Z2)*Z2))) - One,
     ! (le[Z2, X] || le[X, Zero])];
   X = Abs[(((Three/Four) - (Two/Three))*Three) - (One/Four)]; Y = X; Z = X;
   doWhile[Z1 = Z;
     Z = ((One/Two) - (((One/Two) - ((Half*Z1) + ((ThirtyTwo*Z1)*Z1))) + (One/Two))) + (One/Two),
     ! (le[Z1, Z] || le[Z, Zero])];
   doWhile[
     doWhile[Y1 = Y;
       Y = (Half - ((Half - ((Half*Y1) + ((ThirtyTwo*Y1)*Y1))) + Half)) + Half,
       ! (le[Y1, Y] || le[Y, Zero])];
     X1 = X;
     X = (((Half*X1) + ((ThirtyTwo*X1)*X1)) - F9) + F9,
     ! (le[X1, X] || le[X, Zero])];
   If[ne[X1, Y1] || ne[X1, Z1],
     bad[SER, "Disagreements among the values X1, Y1, Z1, respectively " <> fmt[X1] <> ", " <>
        fmt[Y1] <> ", " <> fmt[Z1] <>
        ", are symptoms of inconsistencies introduced by extra-precise evaluation of arithmetic subexpressions."];
     If[eq[X1, U1] || eq[Y1, U1] || eq[Z1, U1], note["That feature is not tested further by this program."]],
     If[ne[Z1, U1] || ne[Z2, U2],
       If[ge[Z1, U1] || ge[Z2, U2],
         bad[FAIL, "Precision test appears to be inconsistent: U1 = " <> fmt[U1] <> ", Z1 - U1 = " <>
            fmt[Z1 - U1] <> "; U2 = " <> fmt[U2] <> ", Z2 - U2 = " <> fmt[Z2 - U2]],
         If[le[Z1, Zero] || le[Z2, Zero],
           note["Because of unusual Radix or exact rational arithmetic a result Z1 = " <> fmt[Z1] <>
              " or Z2 = " <> fmt[Z2] <> " of an extra-precision test appears inconsistent"]];
         If[ne[Z1, Z2] || gt[Z1, Zero],
           X = Z1/U1; Y = Z2/U2; If[gt[Y, X], X = Y]; Q = -Log[X];
           note["Some subexpressions appear to be calculated extra precisely with about " <>
              fmt[Q/Log[Radix]] <> " extra B-digits"]];
         note["That feature is not tested further by this program."]]]];
   pause[]);

(* ------------------------------------------------------------------ *)
(* Milestones 35 - 50: guard digits, rounding, sticky bit            *)
(* ------------------------------------------------------------------ *)

m35[] := (
   ms = 35;
   If[ge[Radix, Two],
     X = W/(Radix*Radix); Y = X + One; Z = Y - X; T = Z + U2; X = T - Z;
     tst[FAIL, eq[X, U2], "Subtraction is not normalized X=Y,X+Z != Y+Z!"];
     If[eq[X, U2], note["Subtraction appears to be normalized, as it should be."]]];
   Y = F9*One; Z = One*F9; X = F9 - Half;
   Y = (Y - Half) - X; Z = (Z - Half) - X;
   X = One + U2; T = X*Radix; R = Radix*X;
   X = T - Radix; X = X - (Radix*U2);
   T = R - Radix; T = T - (Radix*U2);
   X = X*(Radix - One); T = T*(Radix - One);
   If[eq[X, Zero] && eq[Y, Zero] && eq[Z, Zero] && eq[T, Zero], GMult = True,
     GMult = False; tst[SER, False, "* lacks a Guard Digit, so 1*X != X"]];
   Z = Radix*U2; X = One + Z;
   Y = Abs[(X + Z) - (X*X)] - U2;
   X = One - U2;
   Z = Abs[(X - U2) - (X*X)] - U1;
   tst[FAIL, le[Y, Zero] && le[Z, Zero], "* gets too many final digits wrong."];
   Y = One - U2; X = One + U2; Z = One/Y; Y = Z - X;
   X = One/Three; Z = Three/Nine; X = X - Z;
   T = Nine/TwentySeven; Z = Z - T;
   tst[DEF, eq[X, Zero] && eq[Y, Zero] && eq[Z, Zero],
     "Division lacks a Guard Digit, so error can exceed 1 ulp or 1/3 and 3/9 and 9/27 may disagree"];
   Y = F9/One; X = F9 - Half; Y = (Y - Half) - X;
   X = One + U2; T = X/One; X = T - X;
   If[eq[X, Zero] && eq[Y, Zero] && eq[Z, Zero], GDiv = True,
     GDiv = False; tst[SER, False, "Division lacks a Guard Digit, so X/1 != X"]];
   X = One/(One + U2); Y = (X - Half) - Half;
   tst[SER, lt[Y, Zero], "Computed value of 1/1.000..1 >= 1"];
   X = One - U2; Y = One + (Radix*U2);
   Z = X*Radix; T = Y*Radix; R = Z/Radix; StickyBit = T/Radix;
   X = R - X; Y = StickyBit - Y;
   tst[FAIL, eq[X, Zero] && eq[Y, Zero], "* and/or / gets too many last digits wrong"];
   Y = One - U1; X = One - F9; Y = One - Y;
   T = Radix - U2; Z = Radix - BMinusU2; T = Radix - T;
   If[eq[X, U1] && eq[Y, U1] && eq[Z, U2] && eq[T, U2], GAddSub = True,
     GAddSub = False; tst[SER, False, "- lacks Guard Digit, so cancellation is obscured"]];
   If[ne[F9, One] && ge[F9 - One, Zero],
     bad[SER, "comparison alleges (1-U1) < 1 although subtraction yields (1-U1) - 1 = 0, thereby vitiating such precautions against division by zero as ... if (X == 1.0) {.....} else {.../(X-1.0)...}"]];
   If[GMult && GDiv && GAddSub, note["*, /, and - appear to have guard digits, as they should."]]);

m40[] := (
   ms = 40;
   pause[];
   RMult = other; RDiv = other; RAddSub = other;
   RadixD2 = Radix/Two;
   A1 = Two; done = False;
   doWhile[
     AInvrse = Radix;
     doWhile[X = AInvrse; AInvrse = AInvrse/A1, ! ne[fl[AInvrse], AInvrse]];
     done = eq[X, One] || gt[A1, Three];
     If[! done, A1 = Nine + One],
     ! done];
   If[eq[X, One], A1 = Radix];
   AInvrse = One/A1; X = A1; Y = AInvrse; done = False;
   doWhile[
     Z = (X*Y) - Half;
     tst[FAIL, eq[Z, Half], "X * (1/X) differs from 1"];
     done = eq[X, Radix]; X = Radix; Y = One/X,
     ! done];
   Y2 = One + U2; Y1 = One - U2;
   X = OneAndHalf - U2; Y = OneAndHalf + U2;
   Z = (X - U2)*Y2; T = Y*Y1; Z = Z - X; T = T - X;
   X = X*Y2; Y = (Y + U2)*Y1; X = X - OneAndHalf; Y = Y - OneAndHalf;
   If[eq[X, Zero] && eq[Y, Zero] && eq[Z, Zero] && le[T, Zero],
     X = (OneAndHalf + U2)*Y2;
     Y = (OneAndHalf - U2) - U2;
     Z = (OneAndHalf + U2) + U2;
     T = (OneAndHalf - U2)*Y1;
     X = X - (Z + U2);
     StickyBit = Y*Y1;
     S = Z*Y2;
     T = T - Y;
     Y = (U2 - Y) + StickyBit;
     Z = S - ((Z + U2) + U2);
     StickyBit = (Y2 + U2)*Y1;
     Y1 = Y2*Y1;
     StickyBit = StickyBit - Y2;
     Y1 = Y1 - Half;
     Which[
       eq[X, Zero] && eq[Y, Zero] && eq[Z, Zero] && eq[T, Zero] && eq[StickyBit, Zero] && eq[Y1, Half],
         RMult = rounded; note["Multiplication appears to round correctly."],
       eq[X + U2, Zero] && lt[Y, Zero] && eq[Z + U2, Zero] && lt[T, Zero] &&
         eq[StickyBit + U2, Zero] && lt[Y1, Half],
         RMult = chopped; note["Multiplication appears to chop."],
       True, note["* is neither chopped nor correctly rounded."]];
     If[RMult == rounded && ! GMult, note["Multiplication test appears to be inconsistent... PLEASE NOTIFY KARPINKSI!"]],
     note["* is neither chopped nor correctly rounded."]];
   ms = 45;
   Y2 = One + U2; Y1 = One - U2;
   Z = (OneAndHalf + U2) + U2; X = Z/Y2;
   T = (OneAndHalf - U2) - U2; Y = (T - U2)/Y1;
   Z = (Z + U2)/Y2; X = X - OneAndHalf; Y = Y - T; T = T/Y1;
   Z = Z - (OneAndHalf + U2); T = (U2 - OneAndHalf) + T;
   If[! (gt[X, Zero] || gt[Y, Zero] || gt[Z, Zero] || gt[T, Zero]),
     X = OneAndHalf/Y2; Y = OneAndHalf - U2; Z = OneAndHalf + U2;
     X = X - Y; T = OneAndHalf/Y1; Y = Y/Y1; T = T - (Z + U2);
     Y = Y - Z; Z = Z/Y2; Y1 = (Y2 + U2)/Y2; Z = Z - OneAndHalf;
     Y2 = Y1 - Y2; Y1 = (F9 - U1)/F9;
     Which[
       eq[X, Zero] && eq[Y, Zero] && eq[Z, Zero] && eq[T, Zero] && eq[Y2, Zero] &&
         eq[Y1 - Half, F9 - Half],
         RDiv = rounded; note["Division appears to round correctly."];
         If[! GDiv, note["Division test appears to be inconsistent... PLEASE NOTIFY KARPINKSI!"]],
       lt[X, Zero] && lt[Y, Zero] && lt[Z, Zero] && lt[T, Zero] && lt[Y2, Zero] &&
         lt[Y1 - Half, F9 - Half],
         RDiv = chopped; note["Division appears to chop."]]];
   If[RDiv == other, note["/ is neither chopped nor correctly rounded."]];
   BInvrse = One/Radix;
   tst[FAIL, eq[(BInvrse*Radix) - Half, Half], "Radix * ( 1 / Radix ) differs from 1"]);

m50[] := (
   ms = 50;
   tst[FAIL, eq[(F9 + U1) - Half, Half] && eq[(BMinusU2 + U2) - One, Radix - One],
     "Incomplete carry-propagation in Addition"];
   X = One - (U1*U1); Y = One + (U2*(One - U2)); Z = F9 - Half;
   X = (X - Half) - Z; Y = Y - One;
   If[eq[X, Zero] && eq[Y, Zero], RAddSub = chopped; note["Add/Subtract appears to be chopped."]];
   If[GAddSub,
     X = (Half + U2)*U2; Y = (Half - U2)*U2; X = One + X; Y = One + Y;
     X = (One + U2) - X; Y = One - Y;
     If[eq[X, Zero] && eq[Y, Zero],
       X = (Half + U2)*U1; Y = (Half - U2)*U1; X = One - X; Y = One - Y;
       X = F9 - X; Y = One - Y;
       If[eq[X, Zero] && eq[Y, Zero],
         RAddSub = rounded; note["Addition/Subtraction appears to round correctly."],
         note["Addition/Subtraction neither rounds nor chops."]],
       note["Addition/Subtraction neither rounds nor chops."]],
     note["Addition/Subtraction neither rounds nor chops."]];
   S = One;
   X = One + (Half*(One + Half)); Y = (One + U2)*Half;
   Z = X - Y; T = Y - X; StickyBit = Z + T;
   If[ne[StickyBit, Zero], S = Zero; bad[FLW, "(X - Y) + (Y - X) is non zero!"]];
   StickyBit = Zero;
   If[GMult && GDiv && GAddSub && RMult == rounded && RDiv == rounded &&
      RAddSub == rounded && eq[fl[RadixD2], RadixD2],
     X = (Half + U1)*U2; Y = Half*U2; Z = One + Y; T = One + X;
     If[le[Z - One, Zero] && ge[T - One, U2],
       Z = T + Y; Y = Z - X;
       If[ge[Z - T, U2] && eq[Y - T, Zero],
         X = (Half + U1)*U1; Y = Half*U1; Z = One - Y; T = One - X;
         If[eq[Z - One, Zero] && eq[T - F9, Zero],
           Z = (Half - U1)*U1; T = F9 - Z; Q = F9 - Y;
           If[eq[T - F9, Zero] && eq[(F9 - U1) - Q, Zero],
             Z = (One + U2)*OneAndHalf;
             T = ((OneAndHalf + U2) - Z) + U2;
             X = One + (Half/Radix); Y = One + (Radix*U2); Z = X*Y;
             If[eq[T, Zero] && eq[(X + (Radix*U2)) - Z, Zero],
               If[ne[Radix, Two],
                 X = Two + U2; Y = X/Two;
                 If[eq[Y - One, Zero], StickyBit = S],
                 StickyBit = S]]]]]]];
   If[eq[StickyBit, One], note["Sticky bit apparently used correctly."],
     note["Sticky bit used incorrectly or not at all."]];
   tst[FLW, ! (! GMult || ! GDiv || ! GAddSub || RMult == other || RDiv == other || RAddSub == other),
     "lack(s) of guard digits or failure(s) to correctly round or chop (noted above) count as one flaw in the final tally below"]);

(* ------------------------------------------------------------------ *)
(* Milestones 60 - 100                                                *)
(* ------------------------------------------------------------------ *)

m60[] := (
   ms = 60;
   Random9 = Sqrt[Three3]; Random1 = Third; iI = 1;
   doWhile[X = rnd[]; Y = rnd[]; Z9 = Y*X; Z = X*Y; Z9 = Z - Z9; iI = iI + 1,
     ! (iI > NoTrials || ne[Z9, Zero])];
   If[iI == NoTrials,
     Random1 = One + (Half/Three); Random2 = (U2 + U1) + One;
     Z = Random1*Random2; Y = Random2*Random1;
     Z9 = ((One + (Half/Three))*((U2 + U1) + One)) - ((One + (Half/Three))*((U2 + U1) + One))];
   If[! (iI == NoTrials || eq[Z9, Zero]), bad[DEF, "X * Y == Y * X trial fails."],
     note["No failures found in " <> ToString[NoTrials] <> " integer pairs."]]);

m70[] := Module[{stop},
   ms = 70;
   tst[FAIL, eq[Zero, Sqrt[Zero]] && eq[-Zero, Sqrt[-Zero]] && eq[One, Sqrt[One]],
     "Square root of 0.0, -0.0 or 1.0 wrong"];
   MinSqEr = Zero; MaxSqEr = Zero; J = Zero;
   X = Radix; OneUlp = U2; sqXMinX[SER];
   X = BInvrse; OneUlp = BInvrse*U1; sqXMinX[SER];
   X = U1; OneUlp = U1*U1; sqXMinX[SER];
   If[ne[J, Zero], pause[]];
   J = Zero; X = Two; Y = Radix;
   If[ne[Radix, One], doWhile[X = Y; Y = Radix*Y, ! ge[Y - X, NoTrials]]];
   OneUlp = X*U2; iI = 1; stop = False;
   whileG[iI <= NoTrials && ! stop,
     X = X + One; sqXMinX[DEF];
     If[gt[J, Zero], stop = True, iI = iI + 1]];
   (* Test for sqrt monotonicity. *)
   iI = -1; X = BMinusU2; Y = Radix; Z = Radix + (Radix*U2);
   NotMonot = False; Monot = False;
   whileG[! (NotMonot || Monot),
     iI = iI + 1; X = Sqrt[X]; Q = Sqrt[Y]; Z = Sqrt[Z];
     If[gt[X, Q] || gt[Q, Z], NotMonot = True,
       Q = fl[Q + Half];
       If[! (iI > 0 || eq[Radix, Q*Q]), Monot = True,
         If[iI > 0,
           If[iI > 1, Monot = True, Y = Y*BInvrse; X = Y - U1; Z = Y + U1],
           Y = Q; X = Y - U2; Z = Y + U2]]]];
   If[Monot, note["sqrt has passed a test for Monotonicity."],
     bad[DEF, "sqrt(X) is non-monotonic for X near " <> fmt[Y]]]];

m80[] := (
   ms = 80;
   MinSqEr = MinSqEr + Half; MaxSqEr = MaxSqEr - Half;
   Y = (Sqrt[One + U2] - One)/U2;
   SqEr = (Y - One) + (U2/Eight); If[gt[SqEr, MaxSqEr], MaxSqEr = SqEr];
   SqEr = Y + (U2/Eight); If[lt[SqEr, MinSqEr], MinSqEr = SqEr];
   Y = ((Sqrt[F9] - U2) - (One - U2))/U1;
   SqEr = Y + (U1/Eight); If[gt[SqEr, MaxSqEr], MaxSqEr = SqEr];
   SqEr = (Y + One) + (U1/Eight); If[lt[SqEr, MinSqEr], MinSqEr = SqEr];
   OneUlp = U2; X = OneUlp;
   Do[
     Y = Sqrt[((X + U1) + X) + F9];
     Y = ((Y - U2) - ((One - U2) + X))/OneUlp;
     Z = (((((U1 - X) + F9)*Half)*X)*X)/OneUlp;
     SqEr = (Y + Half) + Z; If[lt[SqEr, MinSqEr], MinSqEr = SqEr];
     SqEr = (Y - Half) + Z; If[gt[SqEr, MaxSqEr], MaxSqEr = SqEr];
     If[indx == 1 || indx == 3,
       X = (OneUlp*sgn[X])*fl[Eight/(Nine*Sqrt[OneUlp])],
       OneUlp = U1; X = -OneUlp],
     {indx, 1, 3}]);

m85[] := (
   ms = 85;
   SqRWrng = False; Anomaly = False; RSqrt = other;
   If[ne[Radix, One],
     dD = fl[Half + pw[Radix, (One + prec) - fl[prec]]];
     X = dD/Radix; Y = dD/A1;
     If[ne[X, fl[X]] || ne[Y, fl[Y]], Anomaly = True,
       X = Zero; Z2 = X; Y = One; Y2 = Y; Z1 = Radix - One; FourD = Four*dD;
       doWhile[
         If[gt[Y2, Z2],
           Q = Radix; Y1 = Y;
           doWhile[X1 = Abs[Q + (fl[Half - (Q/Y1)]*Y1)]; Q = Y1; Y1 = X1, ! le[X1, Zero]];
           If[le[Q, One], Z2 = Y2; Z = Y]];
         Y = Y + Two; X = X + Eight; Y2 = Y2 + X;
         If[ge[Y2, FourD], Y2 = Y2 - FourD],
         ! ge[Y, dD]];
       X8 = FourD - Z2; Q = (X8 + (Z*Z))/FourD; X8 = X8/Eight;
       If[ne[Q, fl[Q]], Anomaly = True,
         brk = False;
         doWhile[X = Z1*Z; X = X - (fl[X/Radix]*Radix);
           If[eq[X, One], brk = True, Z1 = Z1 - One],
           ! (brk || le[Z1, Zero])];
         If[le[Z1, Zero] && ! brk, Anomaly = True,
           If[gt[Z1, RadixD2], Z1 = Z1 - Radix];
           doWhile[newD[], ! ge[U2*dD, F9]];
           If[ne[(dD*Radix) - dD, W - dD], Anomaly = True,
             Z2 = dD; iI = 0;
             Y = dD + ((One + Z)*Half); X = (dD + Z) + Q; sr3750[];
             Y = (dD + ((One - Z)*Half)) + dD; X = (dD - Z) + dD; X = (X + Q) + X; sr3750[];
             newD[];
             If[ne[dD - Z2, W - Z2], Anomaly = True,
               Y = (dD - Z2) + (Z2 + ((One - Z)*Half));
               X = (dD - Z2) + ((Z2 - Z) + Q); sr3750[];
               Y = (One + Z)*Half; X = Q; sr3750[];
               If[iI == 0, Anomaly = True]]]]]];
     If[iI == 0 || Anomaly,
       bad[FAIL, "Anomalous arithmetic with Integer < Radix^Precision = " <> fmt[W] <>
          " fails test whether sqrt rounds or chops."];
       SqRWrng = True]];
   If[! Anomaly,
     If[! (lt[MinSqEr, Zero] || gt[MaxSqEr, Zero]),
       RSqrt = rounded; note["Square root appears to be correctly rounded."],
       If[gt[MaxSqEr + U2, U2 - Half] || gt[MinSqEr, Half] || lt[MinSqEr + Radix, Half],
         SqRWrng = True,
         RSqrt = chopped; note["Square root appears to be chopped."]]]];
   If[SqRWrng,
     note["Square root is neither chopped nor correctly rounded. Observed errors run from " <>
        fmt[MinSqEr - Half] <> " to " <> fmt[Half + MaxSqEr] <> " ulps."];
     tst[SER, lt[MaxSqEr - MinSqEr, Radix*Radix], "sqrt gets too many last digits wrong"]]);

m90[] := (
   ms = 90;
   pause[];
   nN = 0; iI = 0; Z = -Zero; mM = 3; brk = False;
   doWhile[
     X = One; sr3980[];
     If[iI <= 10, iI = 1023; sr3980[]];
     If[eq[Z, MinusOne], brk = True, Z = MinusOne; iI = -4],
     ! brk];
   printIfNPositive[];
   N1 = nN; nN = 0; Z = A1;
   mM = Floor[(Two*Log[W])/Log[A1]];
   brk = False;
   doWhile[X = Z; iI = 1; sr3980[];
     If[eq[Z, AInvrse], brk = True, Z = AInvrse],
     ! brk];
   ms = 100;
   mM = NoTrials; Z = Three;
   doWhile[
     X = Z; iI = 1; sr3980[];
     doWhile[Z = Z + Two, eq[Three*fl[Z/Three], Z]],
     lt[Z, Eight*Three]];
   If[nN > 0, note["Errors like this may invalidate financial calculations involving interest rates."]];
   printIfNPositive[];
   nN = nN + N1;
   If[nN == 0, note["... no discrepancies found."]];
   If[nN > 0, pause[]]);

(* ------------------------------------------------------------------ *)
(* Milestones 110 - 220: ADAPTED                                      *)
(* ------------------------------------------------------------------ *)

(* A value that is a positive real number (not Underflow[]/Overflow[]). *)
posQ[x_] := Head[x] === Real && cmpv[x, 0] === 1;
(* exact difference sign of a real against an exact rational (small exponents only) *)

(* Largest k >= 0 with test[Radix^(s k)] True, by doubling then
   bisection; kmax caps the doubling (2^40 >> MPFR's exponent range). *)
bisectK[test_, s_] := Module[{lo = 0, hi = 1, mid},
   While[hi < 2^40 && TrueQ[test[pw[Radix, s*hi]]], lo = hi; hi = 2*hi];
   While[hi - lo > 1, mid = Quotient[lo + hi, 2];
     If[TrueQ[test[pw[Radix, s*mid]]], lo = mid, hi = mid]];
   lo];

(* Evaluate a held expression, recording which messages fired. *)
SetAttributes[withMsgs, HoldFirst];
withMsgs[e_] := Module[{fired = {}, r},
   r = Quiet[Check[e, $chkFailed]];
   If[r === $chkFailed,
     r = Quiet[e];
     (* message names are spelled literally: General::unfl evaluates to
        its text, so it cannot be passed to Check through a variable *)
     fired = Join[
       If[TrueQ[Quiet[Check[e; False, True, {General::unfl}]]], {"General::unfl"}, {}],
       If[TrueQ[Quiet[Check[e; False, True, {General::ovfl}]]], {"General::ovfl"}, {}],
       If[TrueQ[Quiet[Check[e; False, True, {Power::infy}]]], {"Power::infy"}, {}],
       If[TrueQ[Quiet[Check[e; False, True, {Infinity::indet}]]], {"Infinity::indet"}, {}],
       If[TrueQ[Quiet[Check[e; False, True, {Power::indet}]]], {"Power::indet"}, {}]];
     If[fired === {}, fired = {"(a message fired; Check[..., {General::unfl}] does not match it by name inside a package context)"}]];
   {r, fired}];

(* Run PARANOIA's checks from milestone 120 onward on one tier.
   tier: "machine" or "absolute"; ufth = UfThold = E0 (no subnormals);
   vmax = V; kc = exponent used for C. *)
tierChecks[tier_, ufth_, vmax_, cval_] := Module[{tag, x, y, z, r, mf, v9, v0, i, below, back},
   tag = "[" <> tier <> " tier] ";
   UfThold = ufth; E0 = ufth; V = vmax; C0 = cval; CInvrse = One/C0;
   If[lt[Radix, Two], HInvrse = Two, HInvrse = Radix];
   H = One/HInvrse;
   PseudoZero = Zero;
   (* ---- 110/120: the floor itself, and below it ------------------ *)
   ms = 120;
   Z = E0; tstPtUf[];
   {below, mf} = withMsgs[E0*H];
   If[tier === "machine",
     (* Below DBL_MIN a machine result is recomputed as a 53-bit MPFR. *)
     If[posQ[below] && ! MachineNumberQ[below] && cmpv[below, SetPrecision[E0, Infinity]/HInvrse] === 0,
       policy[tag <> "E0 * H = " <> fmt[below] <> " is a positive 53-bit MPFR (no subnormals; values continue as MPFR below DBL_MIN). Underflow is not gradual at this boundary by design."],
       bad[SER, tag <> "E0 * H below the machine range is not the exact positive MPFR value E0/Radix: got " <> fmt[below]]];
     back = below*HInvrse;
     If[! (MachineNumberQ[back] && eq[back, E0]),
       bad[SER, tag <> "(E0 * H) * HInvrse does not return to the machine number E0: got " <> fmt[back]],
       note[tag <> "(E0 * H) * HInvrse demotes back to the machine number E0."]];
     tst[FAIL, lt[below, E0] && gt[below, Zero] && lt[-E0, -below],
       tag <> "comparisons across the machine/MPFR boundary are confused"],
     (* absolute tier *)
     If[Head[below] === Underflow,
       policy[tag <> "E0 * H = Underflow[] (messages " <> fmt[mf] <> "): MPFR has no subnormals; underflow is abrupt at the MPFR exponent floor and yields the symbolic surrogate Underflow[], no trap."],
       If[posQ[below] && lt[below, E0], bad[SER, tag <> "found a positive number " <> fmt[below] <> " below the bisected floor E0"],
         bad[SER, tag <> "E0 * H gave " <> fmt[below] <> ", neither a smaller positive number nor Underflow[]"]]]];
   (* products vs differences: X = E0*(1+U2) != E0; X - E0 must not be 0. *)
   x = E0*(One + U2); z = E0;
   {r, mf} = withMsgs[x - z];
   If[! eq[x, z],
     Which[
       posQ[r] && cmpv[r, If[tier === "machine", SetPrecision[x, Infinity] - SetPrecision[z, Infinity], E0*U2]] === 0,
         note[tag <> "X = E0*(1+U2) != E0 and X - E0 = " <> fmt[r] <> " exactly (differences do not underflow before products)."],
       Head[r] === Underflow,
         If[TrueQ[r == 0],
           bad[FLW, tag <> "X != Z but X - Z == 0 (X - Z is Underflow[] and compares equal to 0)"],
           policy[tag <> "X = E0*(1+U2) != E0 and X - E0 = Underflow[] (messages " <> fmt[mf] <> "): the difference underflows (signalled), it is not a false zero. Differences and products underflow at the same threshold, as with any non-gradual underflow."]],
       eq[r, Zero],
         bad[FLW, tag <> "X = " <> fmt[x] <> " is not equal to Z = " <> fmt[z] <> " yet X - Z yields 0 without signalling Underflow"],
       True,
         bad[SER, tag <> "X - Z near the underflow threshold gave " <> fmt[r] <> " instead of the exact difference"]],
     bad[FAIL, tag <> "E0*(1+U2) compares equal to E0"]];
   (* the C program's "Underflow / UfThold" branch for non-gradual underflow *)
   R = Sqrt[E0/UfThold];
   If[le[R, H], Z = R*UfThold; X = Z*(One + ((R*H)*(One + H))),
     Z = UfThold; X = Z*(One + ((H*H)*(One + H)))];
   {r, mf} = withMsgs[X - Z];
   If[! (eq[X, Z] || ne[r, Zero]),
     bad[FLW, tag <> "X = " <> fmt[X] <> " is not equal to Z = " <> fmt[Z] <> " yet X - Z yields " <> fmt[r]]];
   note[tag <> "The Underflow threshold is " <> fmt[UfThold] <> ", below which calculation may suffer larger Relative error than merely roundoff."];
   Y2 = U1*U1; Y = Y2*Y2; Y2 = Y*U1;
   If[le[Y2, UfThold],
     If[gt[Y, E0], bad[DEF, tag <> "Range is too narrow; U1^5 Underflows."],
       bad[SER, tag <> "Range is too narrow; U1^4 Underflows."]]];
   (* ---- 130: (HInvrse)^(2 log_H UfThold) ------------------------- *)
   ms = 130;
   Y = -fl[Half - ((TwoForty*Log[UfThold])/Log[HInvrse])]/TwoForty;
   Y2 = Y + Y;
   {v9, mf} = withMsgs[Power[HInvrse, Y2]];
   Which[
     Head[v9] === Underflow,
       policy[tag <> "(" <> fmt[HInvrse] <> ")^(" <> fmt[Y2] <> ") = Underflow[] (messages " <> fmt[mf] <> "), no trap: the documented surrogate for a result below MPFR's exponent range."],
     tier === "machine" && posQ[v9] && ! MachineNumberQ[v9],
       If[cmpv[v9, Power[SetPrecision[HInvrse, Infinity], Rationalize[SetPrecision[Y2, Infinity], 0]]] === 0,
         policy[tag <> "(" <> fmt[HInvrse] <> ")^(" <> fmt[Y2] <> ") = " <> fmt[v9] <> ", an exact 53-bit MPFR value far below DBL_MIN rather than an underflow (by design)."],
         bad[DEF, tag <> "(" <> fmt[HInvrse] <> ")^(" <> fmt[Y2] <> ") = " <> fmt[v9] <> " is not the exact power"]],
     ! (ge[v9, Zero] && le[v9, ((Radix + Radix) + E9)*UfThold]),
       bad[SER, tag <> "(" <> fmt[HInvrse] <> ")^(" <> fmt[Y2] <> ") = " <> fmt[v9] <> " is not between 0 and underflow threshold " <> fmt[UfThold]],
     ! gt[v9, UfThold*(One + E9)], note[tag <> "This computed value is O.K."],
     True, bad[DEF, tag <> "(" <> fmt[HInvrse] <> ")^(" <> fmt[Y2] <> ") = " <> fmt[v9] <> " is not between 0 and underflow threshold"]];
   (* ---- 150: Z^Q at four nearly extreme values --------------------- *)
   ms = 150;
   nN = 0; Z = A1; Q = fl[Half - (Log[C0]/Log[A1])]; brk = False;
   doWhile[
     X = CInvrse; Y = pw[Z, Q]; isYeqX[];
     Q = -Q; X = C0; Y = pw[Z, Q]; isYeqX[];
     If[lt[Z, One], brk = True, Z = AInvrse],
     ! brk];
   printIfNPositive[];
   If[nN == 0, note[tag <> "Z^Q at four nearly extreme values: no discrepancies found."]];
   (* ---- 160: overflow ---------------------------------------------- *)
   ms = 160;
   Y = -V;
   {z, mf} = withMsgs[HInvrse*Y];     (* the first product past the threshold *)
   V9 = -Y; V0 = V9;
   If[! eq[-(-Y), Y], bad[FLW, tag <> "-(-Y) differs from Y."]];
   If[tier === "machine",
     If[posQ[-z] && ! MachineNumberQ[z] && cmpv[z, HInvrse*SetPrecision[Y, Infinity]] === 0,
       policy[tag <> "HInvrse * (-V) = " <> fmt[z] <> " is the exact 53-bit MPFR value beyond DBL_MAX (no Infinity, no saturation; values continue as MPFR)."],
       bad[SER, tag <> "overflow past " <> fmt[Y] <> " gives " <> fmt[z] <> " instead of the exact MPFR value"]];
     back = (z*H);
     If[! (MachineNumberQ[back] && eq[back, Y]),
       bad[SER, tag <> "(HInvrse*(-V))*H does not demote back to the machine number -V: got " <> fmt[back]]],
     Which[
       Head[z] === Overflow,
         policy[tag <> "HInvrse * (-V) = Overflow[] (messages " <> fmt[mf] <> "): no trap, no signed saturation value; Overflow[] is the symbolic surrogate (and carries no sign). Saturation value V0 taken as V."],
       posQ[-z] && gt[-z, V],
         bad[SER, tag <> "found a finite number " <> fmt[z] <> " beyond the bisected overflow threshold"],
       True, bad[SER, tag <> "overflow past " <> fmt[Y] <> " shrinks to " <> fmt[z]]]];
   note[tag <> "Overflow threshold is V = " <> fmt[V]];
   {v9, mf} = withMsgs[V*One];
   If[mf =!= {} || ! eq[v9, V], bad[DEF, tag <> "V * 1 = " <> fmt[v9] <> " (messages " <> fmt[mf] <> ")"]];
   {v9, mf} = withMsgs[V/One];
   If[mf =!= {} || ! eq[v9, V], bad[DEF, tag <> "V / 1 = " <> fmt[v9] <> " (messages " <> fmt[mf] <> ")"]];
   (* ---- 170: comparisons -------------------------------------------- *)
   ms = 170;
   If[! (lt[-V, V] && lt[-V0, V0] && lt[-UfThold, V] && lt[UfThold, V]),
     bad[FAIL, tag <> "Comparisons involving +-" <> fmt[V] <> ", +-" <> fmt[V0] <> " and +-" <> fmt[UfThold] <> " are confused by Overflow."]];
   (* ---- 175: sqrt near the underflow threshold ---------------------- *)
   ms = 175;
   Do[
     z = Switch[i, 1, UfThold, 2, E0, 3, PseudoZero];
     If[ne[z, Zero],
       v9 = Sqrt[z]; y = v9*v9;
       If[lt[y/(One - (Radix*E9)), z] || gt[y, (One + (Radix*E9))*z],
         If[gt[v9, U1], bad[SER, ""], bad[DEF, ""]];
         msgs[[-1, 3]] = tag <> "Comparison alleges that what prints as Z = " <> fmt[z] <> " is too far from sqrt(Z) ^ 2 = " <> fmt[y]]],
     {i, 1, 3}];
   (* ---- 180: sqrt near the overflow threshold ----------------------- *)
   ms = 180;
   Do[
     z = If[i == 1, V, V0];
     v9 = Sqrt[z]; x = (One - (Radix*E9))*v9; v9 = v9*x;
     If[lt[v9, (One - ((Two*Radix)*E9))*z] || gt[v9, z],
       If[lt[x, W], bad[SER, ""], bad[DEF, ""]];
       msgs[[-1, 3]] = tag <> "Comparison alleges that Z = " <> fmt[z] <> " is too far from sqrt(Z) ^ 2 (" <> fmt[v9] <> ")"],
     {i, 1, 2}];
   (* ---- 190: range balance ------------------------------------------ *)
   ms = 190;
   pause[];
   x = UfThold*V; y = Radix*Radix;
   If[lt[x*y, One] || gt[x, y],
     If[lt[x*y, U1] || gt[x, y/U1], bad[DEF, tag <> "Badly unbalanced range; UfThold * V = " <> fmt[x] <> " is too far from 1."],
       bad[FLW, tag <> "unbalanced range; UfThold * V = " <> fmt[x] <> " is too far from 1."]]];
   (* ---- 200: X / X ---------------------------------------------------- *)
   ms = 200;
   Do[
     x = Switch[i, 1, F9, 2, One + U2, 3, V, 4, UfThold, 5, Radix];
     y = x;
     v9 = ((y/x) - Half) - Half;
     If[ne[v9, Zero],
       If[eq[v9, -U1] && i < 5, bad[FLW, ""], bad[SER, ""]];
       msgs[[-1, 3]] = tag <> "X / X differs from 1 when X = " <> fmt[x] <> "; instead X / X - 1/2 - 1/2 = " <> fmt[v9]],
     {i, 1, 5}];
   ];

adapted[] := Module[{kMin, kMax, kMachMin, kMachMax, c, v, mf, r, eLo, eHi, t},
   ms = 110;
   (* Absolute tier: exponent bisection on Radix^-k / Radix^k. *)
   kMin = bisectK[posQ, -1];
   kMax = bisectK[posQ, 1];
   note["[absolute tier] smallest positive power Radix^-" <> ToString[kMin] <>
      ", largest finite power Radix^" <> ToString[kMax] <> " (found by exponent bisection)"];
   {r, mf} = withMsgs[Power[Radix, -(kMin + 1)]];
   note["[absolute tier] Radix^-(kMin+1) = " <> fmt[r] <> ", messages " <> fmt[mf]];
   {r, mf} = withMsgs[Power[Radix, kMax + 1]];
   note["[absolute tier] Radix^(kMax+1) = " <> fmt[r] <> ", messages " <> fmt[mf]];
   (* the largest finite number: Radix^kMax * (Radix * (1 - U1)) *)
   v = Power[Radix, kMax]*(Radix*F9);
   If[! posQ[v], v = Power[Radix, kMax]];
   eLo = Power[Radix, -kMin];
   c = Power[Radix, -Min[kMin, kMax]];   (* C with 1/C representable *)
   (* Machine tier: only when the mode's numbers are machine numbers. *)
   If[MachineNumberQ[Radix],
     kMachMin = bisectK[MachineNumberQ[#] && posQ[#] &, -1];
     kMachMax = bisectK[MachineNumberQ[#] && posQ[#] &, 1];
     note["[machine tier] smallest positive machine power Radix^-" <> ToString[kMachMin] <>
        ", largest machine power Radix^" <> ToString[kMachMax]];
     t = Power[Radix, kMachMax]*(Radix*F9);
     If[! MachineNumberQ[t], bad[SER, "[machine tier] Radix^kMachMax * (Radix*(1-U1)) is not a machine number: " <> fmt[t]]];
     UfTholdM = Power[Radix, -kMachMin]; VM = t;
     tierChecks["machine", UfTholdM, VM, UfTholdM],
     note["[machine tier] skipped: Radix " <> fmt[Radix] <> " is not a machine number in this mode"];
     UfTholdM = Missing[]; VM = Missing[]];
   ms = 110;
   tierChecks["absolute", eLo, v, c];
   (* ---- 140: X^((X + 1) / (X - 1)) vs exp(2) (not range related) ------ *)
   m140[];
   (* ---- 210: division by zero --------------------------------------- *)
   ms = 210;
   MyZero = Zero;
   {r, mf} = withMsgs[One/MyZero];
   policy["1 / 0 (One/Zero, Zero = " <> fmt[MyZero] <> ") produces " <> fmt[r] <> ", messages " <> fmt[mf] <> "; no trap."];
   {r, mf} = withMsgs[Zero/MyZero];
   policy["0 / 0 produces " <> fmt[r] <> ", messages " <> fmt[mf] <> "; no trap."];
   {r, mf} = withMsgs[(One - One)/(One - One)];
   policy["(1-1) / (1-1) (a computed zero, " <> fmt[One - One] <> ") produces " <> fmt[r] <> ", messages " <> fmt[mf] <> "."];
   ms = 220;
   <|"UfTholdMachine" -> UfTholdM, "VMachine" -> VM, "UfThold" -> eLo, "V" -> v,
     "kMin" -> kMin, "kMax" -> kMax|>];

m140[] := Module[{iter = 0, go = True},
   ms = 140;
   X = Zero; iI = 2; Y = Two*Three; Q = Zero; nN = 0;
   doWhile[Z = X; iI = iI + 1; Y = Y/(iI + iI); R = Y + Q; X = Z + R; Q = (Z - X) + R, gt[X, Z]];
   Z = (OneAndHalf + (One/Eight)) + (X/(OneAndHalf*ThirtyTwo));
   X = Z*Z; Exp2 = X*X;
   X = F9; Y = X - U1;
   iI = 1;
   While[go && iter < loopCap,
     iter = iter + 1;
     Z = X - BInvrse;
     Z = (X + One)/(Z - (One - BInvrse));
     Q = pw[X, Z] - Exp2;
     If[gt[Abs[Q], TwoForty*U2],
       nN = 1; V9 = (X - BInvrse) - (One - BInvrse);
       bad[DEF, "Calculated " <> fmt[pw[X, Z]] <> " for (1 + (" <> fmt[V9] <> "))^(" <> fmt[Z] <>
          "); differs from correct value " <> fmt[Exp2] <> " by " <> fmt[Q] <>
          ". This much error may spoil financial calculations involving tiny interest rates."];
       go = False,
       Z = ((Y - X)*Two) + Y; X = Y; Y = Z;
       Z = One + ((X - F9)*(X - F9));
       If[gt[Z, One] && iI < NoTrials, iI = iI + 1,
         If[gt[X, One],
           If[nN == 0, note["X^((X + 1) / (X - 1)) vs. exp(2): Accuracy seems adequate."]];
           go = False,
           X = One + U2; Y = U2 + U2; Y = Y + X; iI = 1]]]]];

(* ------------------------------------------------------------------ *)
(* Driver                                                             *)
(* ------------------------------------------------------------------ *)

modeName[m_] := If[StringQ[m], m, "MPFR " <> ToString[m[[2]]]];

Options[ParanoiaRun] = {"Compare" -> "Exact", "Verbose" -> True};

ParanoiaRun[mode_, opts : OptionsPattern[]] := Module[{cnt = 0, res, rng},
   $cmp = OptionValue["Compare"];
   errCnt = {0, 0, 0, 0}; msgs = {}; notes = {}; pol = {}; ms = 0;
   Which[
     mode === "Machine", mk = Function[k, N[k]],
     MatchQ[mode, {"MPFR", _Integer}], With[{d = mode[[2]]}, mk = Function[k, N[k, d]]],
     mode === "Mixed",
       (* even index -> machine, odd -> 40-digit MPFR; One is machine *)
       mk = Function[k, cnt = cnt + 1; If[OddQ[cnt], N[k, 40], N[k]]],
     True, Return[$Failed]];
   makeConstants[];
   m7[]; m10[]; m20[]; m30[]; m35[]; m40[]; m50[]; m60[]; m70[]; m80[]; m85[]; m90[];
   rng = adapted[];
   res = <|"Mode" -> modeName[mode], "Compare" -> $cmp,
      "Radix" -> Radix, "Precision" -> prec,
      "RadixMachineQ" -> MachineNumberQ[Radix], "PrecisionMark" -> Precision[U1],
      "U1" -> U1, "U2" -> U2,
      "RMult" -> RMult, "RDiv" -> RDiv, "RAddSub" -> RAddSub, "RSqrt" -> RSqrt,
      "StickyBit" -> StickyBit, "Range" -> rng,
      "Failure" -> errCnt[[1]], "Serious" -> errCnt[[2]], "Defect" -> errCnt[[3]], "Flaw" -> errCnt[[4]],
      "Messages" -> msgs, "Notes" -> notes, "Policy" -> pol|>;
   If[OptionValue["Verbose"], printReport[res]];
   res];

rname[r_] := Switch[r, 1, "rounded", 2, "chopped", _, "other"];

printReport[res_] := (
   Print["==== PARANOIA: mode ", res["Mode"], ", comparisons ", res["Compare"], " ===="];
   Print["Radix = ", fmt[res["Radix"]], ", Precision = ", fmt[res["Precision"]],
     " (MachineNumberQ: ", res["RadixMachineQ"], ", Precision[U1] = ", res["PrecisionMark"], ")"];
   Print["U1 = ", fmt[res["U1"]], ", U2 = ", fmt[res["U2"]]];
   Print["Multiply ", rname[res["RMult"]], ", divide ", rname[res["RDiv"]], ", add/subtract ",
     rname[res["RAddSub"]], ", sqrt ", rname[res["RSqrt"]], ", sticky bit ", fmt[res["StickyBit"]]];
   Scan[Print["  ", #[[1]], " (milestone ", #[[2]], "): ", #[[3]]] &, res["Messages"]];
   Scan[Print["  POLICY (milestone ", #[[1]], "): ", #[[2]]] &, res["Policy"]];
   Print["FAILUREs ", res["Failure"], ", SERIOUS DEFECTs ", res["Serious"], ", DEFECTs ",
     res["Defect"], ", FLAWs ", res["Flaw"]];
   Print[""]);

(* ------------------------------------------------------------------ *)
(* Gate (used by tests/test_paranoia.c and tests/paranoia/run_paranoia.m) *)
(* ------------------------------------------------------------------ *)

gateModes = {"Machine", {"MPFR", 30}, {"MPFR", 50}, {"MPFR", 100}, "Mixed"};

(* Standing candidate bugs, ratcheted: {mode, class, text prefix}.  A
   FAILURE/SERIOUS that is not listed fails the gate; a listed entry that
   no longer fires is reported as FIXED so the line can be deleted.
   Set the environment variable PARANOIA_STRICT=1 to ignore this list. *)
known = {};
   (* Empty: the first runs' findings -- division by a number in MPFR's
      bottom binades overflowed through its reciprocal, and MPFR division
      rounded twice -- are fixed (v0.6615, one mpfr_div: Divide's MPFR arm
      and the K1 fusion), so the gate is assert-empty. *)

knownQ[mode_, m_] := AnyTrue[known,
   (#[[1]] === All || #[[1]] === mode) && #[[2]] === m[[1]] &&
     StringLength[m[[3]]] >= StringLength[#[[3]]] &&
     StringTake[m[[3]], StringLength[#[[3]]]] === #[[3]] &];

ParanoiaGate[] := Module[{strict, runs, new = {}, seen = {}, fixed, r, mode, gated},
   strict = Environment["PARANOIA_STRICT"] === "1";
   runs = Table[ParanoiaRun[m], {m, gateModes}];
   ParanoiaRun["Machine", "Compare" -> "Language"];   (* informational *)
   Do[
     mode = r["Mode"];
     gated = Select[r["Messages"], MemberQ[{"FAILURE", "SERIOUS"}, #[[1]]] &];
     Do[
       If[! strict && knownQ[mode, m],
         AppendTo[seen, Select[known, (#[[1]] === All || #[[1]] === mode) && #[[2]] === m[[1]] &&
             StringLength[m[[3]]] >= StringLength[#[[3]]] &&
             StringTake[m[[3]], StringLength[#[[3]]]] === #[[3]] &]],
         AppendTo[new, {mode, m[[1]], m[[2]], m[[3]]}]],
       {m, gated}],
     {r, runs}];
   fixed = If[strict, {}, Complement[known, Union[Flatten[seen, 1]]]];
   Print["==== PARANOIA gate ===="];
   Do[Print[r["Mode"], ": FAILURE ", r["Failure"], ", SERIOUS ", r["Serious"],
       ", DEFECT ", r["Defect"], ", FLAW ", r["Flaw"], ", Radix ", fmt[r["Radix"]],
       ", Precision ", fmt[r["Precision"]]], {r, runs}];
   Scan[Print["NEW ", #[[2]], " [", #[[1]], ", milestone ", #[[3]], "]: ", #[[4]]] &, new];
   Scan[Print["FIXED (delete from known): ", fmt[#]] &, fixed];
   Print[If[new === {}, "PARANOIA gate: PASS", "PARANOIA gate: FAIL"]];
   <|"Pass" -> (new === {}), "New" -> new, "Fixed" -> fixed,
     "Counts" -> Table[{r["Mode"], r["Failure"], r["Serious"], r["Defect"], r["Flaw"]}, {r, runs}]|>];

End[];

EndPackage[];
