(* Mathilda — symbolic transformation rules for BesselJ.
 *
 * Numeric evaluation (machine / arbitrary precision, real / complex) lives in
 * the C builtin src/special_functions/besselj.c. The builtin returns NULL for
 * symbolic / exact-numeric arguments, which lets these DownValues fire.
 *
 * 1. Half-integer order -> elementary (spherical Bessel) closed form, but ONLY
 *    for non-numeric (symbolic) argument z. With numeric z the call is left
 *    unevaluated (matching Mathematica: BesselJ[11/2, 1] stays put), so that
 *    machine-precision evaluation goes through the accurate C path rather than
 *    a catastrophically-cancelling elementary expression.
 *
 * 2. Negative integer order reflection: J_{-n}(z) = (-1)^n J_n(z).
 *
 * Every top-level statement is ';'-terminated (newlines do not separate
 * statements in a Get-loaded .m file).
 *)

Unprotect[BesselJ];

(* Upward recurrence for half-integers >= 3/2 (symbolic z only). *)
BesselJ[n_, z_] /; IntegerQ[n - 1/2] && n >= 3/2 && !NumberQ[z] :=
    (2 (n - 1)/z) BesselJ[n - 1, z] - BesselJ[n - 2, z];

(* Downward recurrence for half-integers <= -3/2 (symbolic z only). *)
BesselJ[n_, z_] /; IntegerQ[n - 1/2] && n <= -3/2 && !NumberQ[z] :=
    (2 (n + 1)/z) BesselJ[n + 1, z] - BesselJ[n + 2, z];

(* Base cases. *)
BesselJ[1/2, z_] /; !NumberQ[z] := Sqrt[2/(Pi z)] Sin[z];
BesselJ[-1/2, z_] /; !NumberQ[z] := Sqrt[2/(Pi z)] Cos[z];

(* Negative integer order reflection (fires after the C builtin declines). *)
BesselJ[n_Integer, z_] /; n < 0 := (-1)^n BesselJ[-n, z];

Protect[BesselJ];

(* ------------------------------------------------------------------ *)
(* BesselK -- modified Bessel function of the second kind.            *)
(*                                                                    *)
(* Numeric evaluation (machine / arbitrary precision, real / complex) *)
(* lives in the C builtin bessel.c. As with BesselJ, half-integer     *)
(* closed forms fire only for non-numeric z so inexact arguments take *)
(* the accurate C path (BesselK[11/2, 1] stays unevaluated). Unlike   *)
(* BesselJ, K is EVEN in the order: K_{-nu} = K_nu, and the recurrence *)
(* and half-integer base case carry a '+' instead of a '-'.           *)
(* ------------------------------------------------------------------ *)

Unprotect[BesselK];

(* Upward recurrence for half-integers >= 3/2 (symbolic z only):
   K_{nu+1} = K_{nu-1} + (2 nu / z) K_nu. *)
BesselK[n_, z_] /; IntegerQ[n - 1/2] && n >= 3/2 && !NumberQ[z] :=
    (2 (n - 1)/z) BesselK[n - 1, z] + BesselK[n - 2, z];

(* Even order: reflect negative half-integers to positive (symbolic z). *)
BesselK[n_, z_] /; IntegerQ[n - 1/2] && n <= -1/2 && !NumberQ[z] :=
    BesselK[-n, z];

(* Base case K_{1/2}(z) = Sqrt[Pi/(2 z)] Exp[-z]. *)
BesselK[1/2, z_] /; !NumberQ[z] := Sqrt[Pi/(2 z)] Exp[-z];

(* Even order for integer index: K_{-n} = K_n (fires after the C builtin). *)
BesselK[n_Integer, z_] /; n < 0 := BesselK[-n, z];

Protect[BesselK];

(* ------------------------------------------------------------------ *)
(* BesselI -- modified Bessel function of the first kind.             *)
(*                                                                    *)
(* Numeric evaluation (machine / arbitrary precision, real / complex) *)
(* lives in the C builtin bessel.c. As with BesselJ/BesselK, the       *)
(* half-integer closed forms fire only for non-numeric z so inexact    *)
(* arguments take the accurate C path (BesselI[11/2, 1] stays          *)
(* unevaluated). I is EVEN only for INTEGER order (I_{-n} = I_n); for   *)
(* half-integer order it is NOT (I_{1/2} = Sinh, I_{-1/2} = Cosh), so   *)
(* like BesselJ it needs BOTH base cases and an up/down recurrence.     *)
(* The recurrence (DLMF 10.29.1) is I_{nu-1} - I_{nu+1} = (2 nu/z) I_nu *)
(* (a '-' between the two I terms, vs BesselJ's '+').                   *)
(* ------------------------------------------------------------------ *)

Unprotect[BesselI];

(* Upward recurrence for half-integers >= 3/2 (symbolic z only):
   I_n = I_{n-2} - (2 (n-1)/z) I_{n-1}. *)
BesselI[n_, z_] /; IntegerQ[n - 1/2] && n >= 3/2 && !NumberQ[z] :=
    BesselI[n - 2, z] - (2 (n - 1)/z) BesselI[n - 1, z];

(* Downward recurrence for half-integers <= -3/2 (symbolic z only):
   I_n = I_{n+2} + (2 (n+1)/z) I_{n+1}. *)
BesselI[n_, z_] /; IntegerQ[n - 1/2] && n <= -3/2 && !NumberQ[z] :=
    BesselI[n + 2, z] + (2 (n + 1)/z) BesselI[n + 1, z];

(* Base cases (half-integer order is NOT even, so both are needed). *)
BesselI[1/2, z_] /; !NumberQ[z] := Sqrt[2/(Pi z)] Sinh[z];
BesselI[-1/2, z_] /; !NumberQ[z] := Sqrt[2/(Pi z)] Cosh[z];

(* Even order for integer index: I_{-n} = I_n (fires after the C builtin). *)
BesselI[n_Integer, z_] /; n < 0 := BesselI[-n, z];

Protect[BesselI];

(* ------------------------------------------------------------------ *)
(* BesselY -- Bessel function of the second kind.                     *)
(*                                                                    *)
(* Numeric evaluation (machine / arbitrary precision, real / complex) *)
(* lives in the C builtin bessel.c. As with BesselJ, half-integer     *)
(* closed forms fire only for non-numeric z so inexact arguments take *)
(* the accurate C path (BesselY[11/2, 1] stays unevaluated). Like      *)
(* BesselJ, Y obeys Y_{-n} = (-1)^n Y_n for integer n and the SAME    *)
(* three-term recurrence (a '-' between the two Y terms); only the     *)
(* half-integer base cases differ (Y_{1/2} = -Sqrt[2/(Pi z)] Cos[z],   *)
(* Y_{-1/2} = Sqrt[2/(Pi z)] Sin[z]).                                  *)
(* ------------------------------------------------------------------ *)

Unprotect[BesselY];

(* Upward recurrence for half-integers >= 3/2 (symbolic z only). *)
BesselY[n_, z_] /; IntegerQ[n - 1/2] && n >= 3/2 && !NumberQ[z] :=
    (2 (n - 1)/z) BesselY[n - 1, z] - BesselY[n - 2, z];

(* Downward recurrence for half-integers <= -3/2 (symbolic z only). *)
BesselY[n_, z_] /; IntegerQ[n - 1/2] && n <= -3/2 && !NumberQ[z] :=
    (2 (n + 1)/z) BesselY[n + 1, z] - BesselY[n + 2, z];

(* Base cases (half-integer order is NOT even, so both are needed). *)
BesselY[1/2, z_] /; !NumberQ[z] := -Sqrt[2/(Pi z)] Cos[z];
BesselY[-1/2, z_] /; !NumberQ[z] := Sqrt[2/(Pi z)] Sin[z];

(* Negative integer order reflection (fires after the C builtin declines). *)
BesselY[n_Integer, z_] /; n < 0 := (-1)^n BesselY[-n, z];

Protect[BesselY];
