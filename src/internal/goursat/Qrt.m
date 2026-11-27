(* ========================================================================== *)
(* goursat/Qrt.m -- p = 1/4, 3/4: quartic radicand R = a (t - e)^4 + b with an *)
(* affine order-4 symmetry z -> i z, z = t - e.                                *)
(* goursat.py: quartic_binomial, qrt_reduce.                                   *)
(* ========================================================================== *)

(* If R = a (t - e)^4 + b return {a, b, e}, else None. *)
quarticBinomial[R_, t_] := Module[{c = allCoeffs[R, t], a, e, b},
  If[Length[c] != 5, Return[None, Module]];
  a = c[[1]]; e = -c[[2]]/(4 a);
  b = Expand[R - a (t - e)^4];
  If[FreeQ[b, t] && !isZero[b], {a, b, e}, None]];

(* F dt/y^k on y^4 = a z^4 + b, k in {1, 3}.  F = Sum_j phi_j(z^4) z^j; with     *)
(* x = z^4 the j-th piece is rational when                                       *)
(*   k = 1:  j = 3 (u = y)   -> phi_3 u^2 du/a;   j = 0 (u = y/z) -> -phi_0 u^2 du/(u^4 - a) *)
(*   k = 3:  j = 3 (u = y)   -> phi_3 du/a;       j = 2 (u = y/z) -> -phi_2 du/(u^4 - a)     *)
(* the other two pieces live on lemniscatic curves of genus 1 (the residual).     *)
qrtReduce[F_, R_, t_, y_, k_] := Module[
  {x, u, a, b, e, z, phi, total = 0, log = {}, Fres = 0, jq, rho, ut},
  {a, b, e} = quarticBinomial[R, t];
  z = t - e;
  phi = coords[F, {1, z, z^2, z^3}, z^4, t, x];
  jq = If[k == 1, 0, 2];                          (* the component using u = y/z *)
  Do[If[phi[[j + 1]] =!= 0,
       Which[
         j == 3,
           rho = gK[(phi[[4]] /. x -> (u^4 - b)/a) If[k == 1, u^2, 1]/a];
           ut = y;
           total += ratInt[rho, u] /. u -> ut;
           AppendTo[log, <|"uniformiser" -> ut, "rho" -> rho|>],
         j == jq,
           rho = gK[-(phi[[j + 1]] /. x -> b/(u^4 - a)) If[k == 1, u^2, 1]/(u^4 - a)];
           ut = y/z;
           total += ratInt[rho, u] /. u -> ut;
           AppendTo[log, <|"uniformiser" -> ut, "rho" -> rho|>],
         True,
           Fres += (phi[[j + 1]] /. x -> z^4) z^j]],
     {j, 0, 3}];
  {total, gK[Fres], log}];
