(* ::Package:: *)

(* ========================================================================== *)
(* Divisors.m -- the residue divisors of a normal (reduced) component c Y and  *)
(*               their principal multiples ([B] Prop 2.13 / 2.24; [T] Ch 5    *)
(*               Sec 3, Ch 6 Sec 1).                                          *)
(*                                                                            *)
(* With c = N/E (E normal, Ep = D E the TOWER derivative, Bp the branch-free  *)
(* part of the normal denominator) the places with residue alpha are those    *)
(* over the roots of A = gcd(Bp, N^m G - alpha^m Ep^m) where Y = alpha Ep/N   *)
(* ([T] p.60 with ' = D; Trager's residueIdeal takes Ep as data).  The        *)
(* residues are constants (resultant criterion), so they live in a number     *)
(* field Q(t) (Trager's residueField) and the divisors in the curve over      *)
(* k(t) = fieldFn[{x, lower}, {t, mu, theta}].  With a Q-basis beta_j of the  *)
(* residues (Trager's qBasis), alpha_i = Sum_j b_ij beta_j, the divisors      *)
(*     delta_j = Sum_i b_ij Delta(alpha_i)                                    *)
(* have degree 0 and lie at finite places when every special-place residue    *)
(* vanishes (the case supported here); N_j delta_j = div(u_j) gives the log    *)
(* part Sum_j (beta_j/N_j) Log[u_j] (Prop 2.13).                              *)
(* ========================================================================== *)

brResidueSetup::usage = "brResidueSetup[T, cv, c, nr] gathers the data of the log part of the normal component c Y whose nonzero-residue polynomial is nr[\"P\"] (brNormalResidues): <|\"F0\" -> residue field over Q (fieldQ or fieldNF), \"Fk\" -> the same over k (fieldFn), \"roots\" -> residues (elements of F0), \"basis\" -> Q-basis beta_j, \"coeffs\" -> b, \"N\", \"E\", \"Ep\", \"Bp\" (residue data in theta), \"B\" -> the integral basis descriptor, \"vars\" -> {x, lower...}|>.";
brResidueSetup[T_Association, cv_Association, c_, nr_Association] := Module[
  {th = T["theta"], tc = Together[c], N, E, Ep, En, Bp, t = Unique["Bronstein`Private`brTau"], rf, F0, Fk, qb, vars},
  N = Expand[Numerator[tc]]; E = Expand[Denominator[tc]];
  Ep = Together[brD[T, E]];
  En = Numerator[Together[First[brSplit[T, E]]]];
  Bp = Numerator[Together[En/PolynomialGCD[En, cv["G"]]]];
  rf = residueField[nr["P"], nr["z"], t];
  F0 = rf["F"];
  vars = Prepend[T["lower"], T["x"]];
  Fk = If[F0["Type"] === "NF", fieldFn[vars, {F0["t"], F0["mu"], F0["theta"]}], fieldFn[vars]];
  qb = qBasis[F0, rf["roots"]];
  <|"F0" -> F0, "Fk" -> Fk, "roots" -> rf["roots"], "basis" -> qb["basis"], "coeffs" -> qb["coeffs"],
    "N" -> N, "E" -> E, "Ep" -> Ep, "Bp" -> Expand[Bp], "B" -> cv["B"], "vars" -> vars|>];

brDivisorIdeal::usage = "brDivisorIdeal[rs, b] is the ideal (over rs[\"Fk\"], on the basis rs[\"B\"]) of the divisor Sum_i b_i Delta(alpha_i), alpha_i = rs[\"roots\"] (Trager's divisorIdeal with the tower-derivative residue data of rs).";
brDivisorIdeal[rs_Association, b_List] :=
  divisorIdeal[rs["Fk"], rs["B"], rs["roots"], b, rs["N"], rs["Ep"], rs["Bp"]];

brPlaceDegree::usage = "brPlaceDegree[rs, alpha] is the degree over k(t) of Delta(alpha): the theta-degree of A = gcd(Bp, N^m G - alpha^m Ep^m).";
brPlaceDegree[rs_Association, alpha_] := Module[{F = rs["Fk"], B = rs["B"], x = rs["B"]["x"]},
  kDeg[F, kGCD[F, rs["Bp"], kc[F, rs["N"]^B["m"] B["G"] - alpha^B["m"] rs["Ep"]^B["m"]], x], x]];

brPrincipalGenerator::usage = "brPrincipalGenerator[rs, b, n] is a generator u (coordinates on rs[\"B\"] over rs[\"Fk\"]) of n delta, delta = Sum_i b_i Delta(alpha_i), or None when n delta is not principal (Trager's exact test, [T] Ch 6 Sec 1: a basis normal at infinity).";
brPrincipalGenerator[rs_Association, b_List, n_Integer] :=
  principalGenerator[rs["Fk"], rs["B"], brDivisorIdeal[rs, n b]];

(* -------------------------------------------------------------------------- *)
(* Divisors with places above the special points (genus 0).                  *)
(*                                                                            *)
(* The residues at the places above theta = infinity (logarithm, Def 2.10)    *)
(* or above 0 and infinity (exponential, Def 2.21) join the finite ones; the  *)
(* Q-basis is taken over all of them.  The divisors are built on the shifted  *)
(* curve theta = a + 1/u (a an ordinary rational point, a != 0): Yh = Y u^k,  *)
(* k = Ceiling[g/m], Yh^m = Gu(u) = u^(k m) G(a + 1/u), so every special place *)
(* is FINITE (above u = 0 for infinity, u = -1/a for theta = 0) and no place  *)
(* of the divisor lies above u = infinity (theta = a is ordinary): Trager's   *)
(* principality test applies.  Finite residue places keep their description  *)
(* Yh = alpha Eh/Nh with Eh/Nh = u^k Ep(a + 1/u)/N(a + 1/u) (residueIdeal on  *)
(* the frame).  An unramified special place is the point (u - u0, Yh - w);    *)
(* a totally ramified one Trager's placeIdealAt.  Only genus 0 is supported:  *)
(* there every degree-0 divisor is principal (the torsion order is 1); for     *)
(* higher genus the specialized torsion computation would have to carry the    *)
(* special places too.                                                        *)
(* -------------------------------------------------------------------------- *)

brSpecialLocations::usage = "brSpecialLocations[T, cv, spec] lists the special points carrying places of nonzero residue, from the special-residue data spec (brLogInfinityResidues / brExpSpecialResidues): {<|\"where\" -> \"inf\" | \"zero\", \"ramified\" -> bool, \"kappa\", \"u\", \"coef\"|>, ...} (coef: lambda for a logarithm, c_w for an exponential).";
brSpecialLocations[T_Association, cv_Association, spec_Association] := If[T["kind"] === "Log",
  If[TrueQ[spec["ramified"]] || spec["lambda"] === 0, {},
    {<|"where" -> "inf", "ramified" -> False, "kappa" -> spec["kappa"], "u" -> spec["u"], "coef" -> spec["lambda"]|>}],
  Join[
    {<|"where" -> "zero", "ramified" -> TrueQ[spec["zero"]["ramified"]], "kappa" -> Lookup[spec["zero"], "kappa", 1],
       "u" -> Lookup[spec["zero"], "u", 1], "coef" -> Lookup[spec["zero"], "c", 0]|>},
    {<|"where" -> "inf", "ramified" -> TrueQ[spec["inf"]["ramified"]], "kappa" -> Lookup[spec["inf"], "kappa", 1],
       "u" -> Lookup[spec["inf"], "u", 1], "coef" -> Lookup[spec["inf"], "c", 0]|>}]];

brPointPlaceIdeal::usage = "brPointPlaceIdeal[F, B, A, w] is the prime ideal of the unramified place (A = 0, Y = w) of the curve of B, A a polynomial of degree 1 in the curve variable and w a scalar of F with w^m = G at the root of A (checked): spanned by A w_0 and w_k - (w^k/d_k mod A) w_0.";
brPointPlaceIdeal[F_Association, B_Association, A_, w_] := Module[{x = B["x"], m = B["m"], rows, dinv},
  If[!kZeroQ[F, kRem[F, kc[F, w^m - B["G"]], A, x]], throwInternal["brPointPlaceNotOnCurve", <|"w" -> w|>]];
  rows = Prepend[
    Table[
      dinv = kExtGCD[F, kRem[F, B["d"][[j + 1]], A, x], A, x];
      If[dinv[[1]] =!= 1, throwInternal["brPointPlaceRamified"]];
      ReplacePart[Table[0, {m}], {1 -> kRem[F, -kc[F, w^j] dinv[[2]], A, x], j + 1 -> 1}],
      {j, 1, m - 1}],
    ReplacePart[Table[0, {m}], 1 -> kc[F, A]]];
  <|"M" -> hnf[F, rows, x], "den" -> 1|>];

brShiftFrame::usage = "brShiftFrame[T, cv, rs, a, u] is the frame theta = a + 1/u of the component: <|\"u\" -> u, \"a\" -> a, \"k\" -> k, \"Bu\" -> basis of Yh^m = Gu, \"Nh\", \"Eh\", \"Bph\" (finite residue data in u)|>.";
brShiftFrame[T_Association, cv_Association, rs_Association, a_, u_Symbol] := Module[
  {th = T["theta"], m = cv["m"], G = cv["G"], g, k, Gu, r, Bph},
  g = Exponent[G, th]; k = Ceiling[g/m];
  Gu = Expand[Together[u^(k m) (G /. th -> a + 1/u)]];
  r = Together[u^k (rs["Ep"] /. th -> a + 1/u)/(rs["N"] /. th -> a + 1/u)];
  Bph = If[FreeQ[rs["Bp"], th], 1, Expand[Numerator[Together[rs["Bp"] /. th -> a + 1/u]]]];
  <|"u" -> u, "a" -> a, "k" -> k, "Bu" -> curveBasis[m, Gu, u],
    "Nh" -> Expand[Denominator[r]], "Eh" -> Expand[Numerator[r]], "Bph" -> Bph|>];

brSpecialPlaceIdeals::usage = "brSpecialPlaceIdeals[F, fr, loc, rhos] is the list {{ideal, rho}, ...} of the place ideals on the frame fr at the special location loc: one per root rho (of z^m - kappa) for an unramified location (Yh = rho u0 scaled to the frame), or the single totally ramified place {ideal, None}.";
brSpecialPlaceIdeals[F_Association, fr_Association, loc_Association, rhos_List] := Module[
  {u = fr["u"], a = fr["a"], k = fr["k"], B = fr["Bu"], pt},
  pt = If[loc["where"] === "inf", 0, -1/a];
  If[loc["ramified"],
    If[GCD[rootMultiplicity[F, B["G"], pt, u], B["m"]] =!= 1,
      throwDecline["several ramified places above a special point"]];
    {{placeIdealAt[F, B, pt], None}},
    {brPointPlaceIdeal[F, B, kc[F, u - pt], kc[F, # loc["u"] If[loc["where"] === "inf", 1, (-1/a)^k]]], #} & /@ rhos]];

brSpecialTaus::usage = "brSpecialTaus[T, cv, locs, rhoLists] is the list of residues tau (one per place, in the order of brSpecialPlaceIdeals) at the special places: logarithm, infinity: -rho lambda; exponential ([B] Def 2.21, pbar the first place above 0): eps r (c_p - c_pbar) with c_p = rho c_w at an unramified place and 0 at a ramified one.";
brSpecialTaus[T_Association, cv_Association, locs_List, rhoLists_List] := Module[{m = cv["m"], cbar, val, out = {}},
  If[T["kind"] === "Log",
    Return[Flatten[MapThread[Function[{loc, rhos}, (-# loc["coef"]) & /@ rhos], {locs, rhoLists}]]]];
  (* exponential: values c_p, then tau = eps r (c_p - c_pbar) *)
  val[loc_, rho_] := If[loc["ramified"], 0, rho loc["coef"]];
  cbar = With[{z = First[locs]}, val[z, If[z["ramified"], None, First[rhoLists[[1]]]]]];
  Do[
    With[{loc = locs[[i]], eps = If[locs[[i]]["where"] === "zero", 1, -1]},
      With[{r = If[loc["ramified"], m/GCD[m, If[loc["where"] === "zero",
                    brPlaceDataZero[T, cv]["e"], brPlaceDataInf[T, cv]["g"]]], 1]},
        If[loc["ramified"], AppendTo[out, eps r (0 - cbar)],
          Do[AppendTo[out, eps r (val[loc, rho] - cbar)], {rho, rhoLists[[i]]}]]]],
    {i, Length[locs]}];
  out];
