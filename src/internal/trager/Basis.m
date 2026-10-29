(* ::Package:: *)

(* ========================================================================== *)
(* Basis.m -- the integral basis of Q(x)[Y]/(Y^m - G), [T] Ch 2 Sec 5 p.30;  *)
(*            [B] Sec 2.2 eq. (11).                                           *)
(*                                                                            *)
(* With G = c Prod p_j^e_j, 0 < e_j < m:                                      *)
(*     w_i = Y^i / d_i,   d_i = Prod_j p_j^Floor[i e_j / m],   i = 0..m-1.     *)
(* The basis is diagonal for the derivation,                                  *)
(*     w_i' = s_i w_i,   s_i = i G'/(m G) - d_i'/d_i,                         *)
(* and multiplicatively closed up to polynomial factors,                      *)
(*     w_i w_j = G^q d_k/(d_i d_j) w_k,  k = (i+j) mod m, q = (i+j) div m,    *)
(* the factor being a polynomial ([T] p.31).  It is normal at infinity        *)
(* ([T] p.30, [B] p.17): with M = deg G the elements                          *)
(*     v_i = x^(-delta_i) w_i,   delta_i = Ceiling[i M/m] - deg d_i,          *)
(* form a local integral basis at infinity.  Because d_1 = 1, the generator   *)
(* itself is w_1 = Y.                                                         *)
(* ========================================================================== *)

curveBasis::usage = "curveBasis[m, G, x] is the basis descriptor of the curve Y^m = G: an Association with \"m\", \"G\", \"x\", \"deg\" (= deg G), \"factors\" ({{p_j, e_j}}), \"d\" ({d_0..d_(m-1)}), \"s\" ({s_0..s_(m-1)}, the logarithmic derivatives of the w_i), \"delta\" ({delta_i}, the exponents of the local integral basis at infinity) and \"mult\" (an m x m table of {k, factor} with w_i w_j = factor w_k).";
curveBasis[m_Integer, G_, x_Symbol] := Module[{pairs, d, s, M, delta, mult},
  pairs = Select[FactorList[G], !FreeQ[#[[1]], x] &];
  d = Table[Expand[Times @@ (#[[1]]^Floor[i #[[2]]/m] & /@ pairs)], {i, 0, m - 1}];
  s = Table[Together[i D[G, x]/(m G) - D[d[[i + 1]], x]/d[[i + 1]]], {i, 0, m - 1}];
  M = Exponent[G, x];
  delta = Table[Ceiling[i M/m] - Exponent[d[[i + 1]], x], {i, 0, m - 1}];
  mult = Table[
    With[{k = Mod[i + j, m], q = Quotient[i + j, m]},
      {k, Together[G^q d[[k + 1]]/(d[[i + 1]] d[[j + 1]])]}],
    {i, 0, m - 1}, {j, 0, m - 1}];
  <|"m" -> m, "G" -> Expand[G], "x" -> x, "deg" -> M, "factors" -> pairs,
    "d" -> d, "s" -> s, "delta" -> delta, "mult" -> mult|>];
