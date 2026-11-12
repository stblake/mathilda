(* Elliptic integrals ON the branch cut 1 - m Sin[phi]^2 = -1 (exactly).
   Found by the integration shootout (2026-11-12): the cubic Legendre form of
   Integrate[1/((x - 1) Sqrt[1 + x^3]), x] with the old root ordering reached
   this point at the verification sample x = 1/2. *)
ph = ArcSin[Sqrt[1 + I Sqrt[3]]]; m = (1 - I Sqrt[3])/2;
Print["1 - m Sin[phi]^2 = ", Simplify[1 - m Sin[ph]^2]];                 (* -1 *)
Print["F machine   ", EllipticF[N[ph], N[m]]];                           (* 1.86639 + 0.790708 I  correct *)
Print["F 30 digits ", N[EllipticF[ph, m], 30]];                          (* WRONG: 0.39373 + 1.46941 I *)
Print["Pi numeric  ", EllipticPi[N[-I/Sqrt[3], 30], N[ph, 30], N[m, 30]]]; (* 2.31183 - 0.240248 I  correct *)
Print["Pi exact    ", TimeConstrained[N[EllipticPi[-I/Sqrt[3], ph, m], 20], 10]];   (* answers *)
(* the same phi, spelled as the substitution x -> 1/2 left it *)
ph2 = ArcSin[Sqrt[(-1 + (-1)^(2/3))/(1/2 + (-1)^(2/3))]];
Print["same phi?   ", PossibleZeroQ[ph2 - ph]];
Print["Pi exact 2  ", TimeConstrained[N[EllipticPi[-I/Sqrt[3], ph2, m], 20], 10]];  (* HANGS: $Aborted *)
