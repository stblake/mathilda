(* Gaussian exponent offsets (campaign step 4 over k(i)): E^(I (k + u)) beside
   E^(I u) is E^(k I) (E^(I u)): the decision writes E^(a + b I) as E^a Z^b with
   Z a symbol for E^I (Lindemann-Weierstrass: E, E^I independent).  F mixes
   shifted exponentials, shifted tangents, literal E^I-type constants and
   constants spelled with u (E^(I u - I (1 + u)) = E^-I, the spelling that made
   a false "no" while the rule was written).  D[F] is elementary, so a False is
   a false certificate.  Second half: random integrands, decisions counted, a
   False contradicted by a verified antiderivative reported. *)
ek = {E^(I u), E^(I (1 + u)), E^(-I (1 + u)), E^(I (2 + u)), E^(2 I u), E^(I + 2 I u), E^u, E^(1 + u)};
tk = {Tan[1 + u], Tan[2 + u], Tan[u], Tanh[1 + u]};
ck = {1, E^I, E^(-I), E^(I u - I (1 + u)), E, 2 I, 1 + I};
mk := {RandomChoice[ek], RandomChoice[tk], Log[u + I], Log[RandomChoice[tk] + I], u, 1, 1/(1 + RandomChoice[ek]), 1/(RandomChoice[tk] - I)};
Do[SeedRandom[sd];
 res = Table[Module[{F},
    F = Sum[RandomChoice[ck] u^RandomInteger[{0, 1}] RandomChoice[mk] RandomChoice[mk]/RandomChoice[{1, u, u + I, 1 + RandomChoice[ek]}], {RandomInteger[{1, 2}]}];
    {F, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[D[F, u], u]], 30, $TO]}], {80}];
 Print["seed ", sd, ": True=", Count[res, {_, True}], " False=", Count[res, {_, False}], " TO=", Count[res, {_, $TO}], " undec=", Count[res, {_, _Risch`ElementaryIntegralQ}]];
 Scan[Print["  FALSE CERT: ", InputForm[#[[1]]]] &, Select[res, #[[2]] === False &]], {sd, {1701, 1702}}];
SeedRandom[1710];
res2 = Table[Module[{f},
    f = RandomChoice[ck] RandomChoice[{1, u, 1/u, 1/(u + I), RandomChoice[tk]}] RandomChoice[ek] RandomChoice[{1, RandomChoice[tk], RandomChoice[ek], u}];
    {f, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[f, u]], 30, $TO], TimeConstrained[Quiet[Integrate`RischTranscendental[f, u]], 30, $TO]}], {80}];
Print["random: True=", Count[res2, {_, True, _}], " False=", Count[res2, {_, False, _}], " TO=", Count[res2, {_, $TO, _}], " undec=", Count[res2, {_, _Risch`ElementaryIntegralQ, _}]];
(* only an ELEMENTARY antiderivative contradicts a False: an Ei/Erf answer is
   the non-elementary form the False asserts *)
sfQ[e_] := !FreeQ[e, ExpIntegralEi | ExpIntegralE | Erf | Erfi | Erfc | LogIntegral | SinIntegral | CosIntegral | Gamma | PolyLog | Hypergeometric2F1];
Scan[If[#[[2]] === False && Head[#[[3]]] =!= Integrate`RischTranscendental && #[[3]] =!= $TO && !sfQ[#[[3]]] && PossibleZeroQ[D[#[[3]], u] - #[[1]]], Print["  CONTRADICTION: ", InputForm[#[[1]]]]] &, res2];
