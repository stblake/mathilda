(* D, Dt and Derivative are now implemented in C (src/deriv.c) and
   registered during core_init(). The rule-based bootstrap that used
   to live in src/internal/deriv.m is no longer loaded. *)

(* The CRC integral table (CRCMathTablesIntegrals.m) is loaded lazily
   the first time Integrate's CRCTable stage is invoked — see
   src/integrate.c's try_crctable().  We do not Get it here so startup
   cost stays low for users who never call Integrate. *)

(* Symbolic transformation rules for BesselJ (half-integer -> elementary,
   negative-integer-order reflection). Numeric evaluation is in the C
   builtin; these rules only fire for symbolic / exact arguments. *)
Get["./src/internal/bessel.m"];
