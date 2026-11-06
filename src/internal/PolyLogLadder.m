(* PolyLogLadder.m -- polylogarithmic antiderivatives of R(x) Log[a x + b]^k.

   Loaded lazily by the Integrate cascade stage "PolyLogLadder"
   (src/calculus/integrate_polylog.c).  All symbols are fully qualified in the
   PolyLogLadder` context; the file is not a package.

   CLASS.  f = R0(x) + Sum_i Sum_k R_ik(x) Log[w_i]^k, where every R is a
   rational function of x and every w_i = a_i x + b_i is LINEAR in x (a single
   log kernel per term: products of two distinct log kernels decline).  The
   denominators may have algebraic or complex roots of any multiplicity (degree
   <= 4, so the roots are radicals), and the parameters may be symbolic.

   METHOD.  Work in the kernel variable W = w_i itself (x = (W - b)/a): then
   Log[w_i] IS Log[W] -- no splitting of a logarithm into a sum, so no branch
   bookkeeping -- and the antiderivative of R(W) Log[W]^k is assembled from
   partial fractions over the roots of R's denominator by four exact identities:

     INT W^n Log[W]^k  = Sum_j (-1)^j k!/(k-j)! W^(n+1) Log[W]^(k-j)/(n+1)^(j+1)
     INT Log[W]^k / W  = Log[W]^(k+1)/(k+1)
     INT Log[W]^k/(W - c) = Log[W]^k Log[1 - W/c]
                            + k Sum_t (-1)^t (k-1)!/(k-1-t)! Log[W]^(k-1-t)
                                    PolyLog[t+2, W/c]                 (c != 0)
     INT Log[W]^k (W - c)^-m = Log[W]^k (W-c)^(1-m)/(1-m)
                               - k/(1-m) INT Log[W]^(k-1) (W-c)^(1-m)/W   (m >= 2)

   (the last recurses on a lower log power).  The third is the polylogarithm
   ladder: d/dW PolyLog[s+1, W/c] = PolyLog[s, W/c]/W and
   d/dW Log[1 - W/c] = 1/(W - c), so every term differentiates back exactly.
   The C stage re-checks the result by differentiation regardless.

   TRIG ROUTE.  A polynomial in x times a rational function of trig/hyperbolic
   functions of linear arguments (x^2 Csc[x], x/(1 + E^x) Cos[x], ...) is
   rewritten with TrigToExp and handed to the exponential-tower Risch stage,
   whose polylogarithm ladder (cherry_polylog_exp.c) closes it.  Without the
   rewrite that ladder never sees trig input.

   Returns the antiderivative, or $Failed when the integrand is outside the
   class (the stage then declines and the cascade continues). *)

PolyLogLadder`$maxDenDegree = 4;

(* INT W^n Log[W]^k dW, integer n != -1. *)
PolyLogLadder`PowLog[n_, k_, W_] :=
  Sum[(-1)^j k!/(k - j)! W^(n + 1) Log[W]^(k - j)/(n + 1)^(j + 1), {j, 0, k}];

(* INT Log[W]^k/(W - c) dW, c != 0: the polylogarithm ladder. *)
PolyLogLadder`SimpleLog[c_, k_, W_] :=
  Log[W]^k Log[1 - W/c] +
    k Sum[(-1)^t (k - 1)!/(k - 1 - t)! Log[W]^(k - 1 - t) PolyLog[t + 2, W/c], {t, 0, k - 1}];

PolyLogLadder`ZeroQ[e_] := TrueQ[Simplify[e] == 0];

(* Distinct roots of the polynomial d in W with multiplicities, or $Failed.
   Factor first: each irreducible factor of degree <= 4 is solved in radicals and
   contributes its roots with the factor's exponent as multiplicity. *)
PolyLogLadder`RootsMult[d_, W_] := Catch[Module[{fl, out = {}, sol, rs},
  If[! PolynomialQ[d, W] || Exponent[d, W] < 1, Return[$Failed]];
  fl = FactorList[d];
  Do[
    With[{p = fe[[1]], e = fe[[2]]},
      If[FreeQ[p, W], Continue[]];
      If[Exponent[p, W] > PolyLogLadder`$maxDenDegree, Throw[$Failed, PolyLogLadder`fail]];
      sol = Quiet[Solve[p == 0, W]];
      If[! ListQ[sol] || Length[sol] != Exponent[p, W], Throw[$Failed, PolyLogLadder`fail]];
      rs = W /. sol;
      If[! FreeQ[rs, Root | W], Throw[$Failed, PolyLogLadder`fail]];
      Do[AppendTo[out, {r, e}], {r, rs}]],
    {fe, fl}];
  If[Total[out[[All, 2]]] =!= Exponent[d, W], Return[$Failed]];
  out], PolyLogLadder`fail];

(* INT R(W) Log[W]^k dW for rational R. *)
PolyLogLadder`Ladder[R_, W_, k_] := Catch[Module[
  {rr, num, den, poly, rem, roots, res = 0, q, h, cb, sub},
  If[k == 0,
    sub = Integrate[R, W];
    Return[If[FreeQ[sub, Integrate], sub, $Failed]]];
  rr = Together[R];
  num = Numerator[rr]; den = Denominator[rr];
  If[! PolynomialQ[num, W] || ! PolynomialQ[den, W], Return[$Failed]];
  poly = PolynomialQuotient[num, den, W];
  rem = PolynomialRemainder[num, den, W];
  (* polynomial part *)
  poly = Expand[poly];
  Do[
    If[Coefficient[poly, W, n] =!= 0,
      res += Coefficient[poly, W, n] PolyLogLadder`PowLog[n, k, W]],
    {n, 0, Exponent[poly, W]}];
  If[PolyLogLadder`ZeroQ[rem], Return[res]];
  If[Exponent[den, W] == 0, Return[res + PolyLogLadder`PowLog[0, k, W] rem/den]];
  roots = PolyLogLadder`RootsMult[den, W];
  If[roots === $Failed, Return[$Failed]];
  Do[
    With[{b = rt[[1]], m = rt[[2]]},
      q = PolynomialQuotient[den, Expand[(W - b)^m], W];
      If[! PolyLogLadder`ZeroQ[PolynomialRemainder[den, Expand[(W - b)^m], W]], Throw[$Failed, PolyLogLadder`fail]];
      h = rem/q;
      Do[
        cb = Simplify[(D[h, {W, m - j}] /. W -> b)/(m - j)!];
        If[! PolyLogLadder`ZeroQ[cb],
          res += cb PolyLogLadder`PoleTerm[b, j, k, W];
          If[! FreeQ[res, $Failed], Throw[$Failed, PolyLogLadder`fail]]],
        {j, 1, m}]],
    {rt, roots}];
  res], PolyLogLadder`fail];

(* INT Log[W]^k/(W - b)^j dW. *)
PolyLogLadder`PoleTerm[b_, j_, k_, W_] := Module[{sub},
  Which[
    PolyLogLadder`ZeroQ[b] && j == 1, Log[W]^(k + 1)/(k + 1),
    PolyLogLadder`ZeroQ[b], PolyLogLadder`PowLog[-j, k, W],
    j == 1, PolyLogLadder`SimpleLog[b, k, W],
    True,
      sub = PolyLogLadder`Ladder[(W - b)^(1 - j)/W, W, k - 1];
      If[sub === $Failed, $Failed,
        Log[W]^k (W - b)^(1 - j)/(1 - j) - k/(1 - j) sub]]];

(* LOG ARGUMENTS BEYOND LINEAR.

   Splitting.  Log[P], P a polynomial of degree >= 2 with real numeric
   coefficients and no real root, is a sum of linear logs over its (non-real,
   conjugate-paired) roots r.  For real x and any real point p,
     Log[P(x)] = Log[P(p)] + Sum_r m_r Log[(x - r)/(p - r)]
   exactly: conjugate factors have arguments summing to exactly 0, and
   P(x)/P(p) > 0.  So Log[1 + x^2] = Log[1 - I x] + Log[1 + I x] (p = 0).  The
   point p is chosen per term of R's partial fractions: a term with its pole at
   a rational p is split at that p, so every kernel equals 1 at the pole and the
   ladder's Log[W] Log[1 - W] is 0 there instead of two divergent logs whose
   coefficients cancel only after the values are simplified (Int_0^1
   Log[1 + x + x^2]/x = Pi^2/9 depends on it); every other term is split at
   p = 0 (Log[P(0)] = Log of the constant term).

   Moebius kernels.  A kernel Log[M], M = (al x + be)/(ga x + de), is integrated
   in u = M: x = (de u - be)/(al - ga u) is Moebius again, so R(x) dx stays a
   rational R(x(u)) x'(u) du and the ladder applies in u.  Replacing (M, Log[M])
   by (1/M, -Log[M]) is exact too (d/dx -Log[M] = (1/M)'/(1/M)); the orientation
   with |u| <= 1 at the stage's sample point 37/100 is taken, which keeps
   PolyLog[n, u/c] off its cut on [0, 1] for the common kernels.

   Polarisation.  A product of two distinct linear kernels is
     L_i L_j = (L_i^2 + L_j^2 - (L_i - L_j)^2)/2,
   and L_i - L_j = Log[w_i] - Log[w_j] differs from Log[w_i/w_j] by a constant
   (2 Pi I n, piecewise), so it is handled as the Moebius kernel u = w_i/w_j
   with Log[u] replaced by L_i - L_j in the ladder's answer: the ladder uses
   only d/du Log[u] = 1/u, which L_i - L_j satisfies exactly, so the answer
   differentiates back to R (L_i - L_j)^2 everywhere, with no branch to track.
   Int_0^1 Log[1 - x] Log[1 + x]/x = -5 Zeta[3]/8 this way. *)

(* {roots with multiplicities} of a splittable P, or $Failed *)
PolyLogLadder`SplitRoots[u_, x_] := Module[{p = Expand[u], cl, rts},
  If[! PolynomialQ[p, x] || Exponent[p, x] < 2, Return[$Failed]];
  cl = CoefficientList[p, x];
  If[! AllTrue[cl, NumericQ[#] && TrueQ[Im[#] == 0] &], Return[$Failed]];
  rts = PolyLogLadder`RootsMult[p, x];
  If[rts === $Failed, Return[$Failed]];
  If[! AllTrue[rts, TrueQ[Abs[Im[N[#[[1]], 30]]] > 10^-20] &], Return[$Failed]];
  rts];

(* Log[P(x)] split at the real point p0 *)
PolyLogLadder`SplitAt[u_, rts_, x_, p0_] :=
  Log[u /. x -> p0] + Total[(#[[2]] Log[Together[(x - #[[1]])/(p0 - #[[1]])]]) & /@ rts];

(* the rational real pole of a partial-fraction term, or None *)
PolyLogLadder`TermPole[t_, x_] := Module[{d = Denominator[Together[t]], fl},
  If[FreeQ[d, x], Return[None]];
  fl = Select[FactorList[d], ! FreeQ[#[[1]], x] &];
  If[Length[fl] == 1 && Exponent[fl[[1, 1]], x] == 1,
    -Coefficient[fl[[1, 1]], x, 0]/Coefficient[fl[[1, 1]], x, 1], None]];

(* f with every splittable Log[P] replaced by linear logs, split per pole *)
PolyLogLadder`SplitAll[f_, x_] := Catch[Module[{logs, sp, lam, g, others, osym, out = 0, terms},
  logs = Union[Cases[f, Log[u_] /; ! FreeQ[u, x], {0, Infinity}]];
  sp = Select[{#, PolyLogLadder`SplitRoots[#[[1]], x]} & /@
         Select[logs, PolynomialQ[#[[1]], x] && Exponent[#[[1]], x] >= 2 &], #[[2]] =!= $Failed &];
  If[sp === {}, Throw[f, PolyLogLadder`split]];
  lam = Table[Unique["PolyLogLadder`S"], {Length[sp]}];
  others = Complement[logs, sp[[All, 1]]];
  osym = Table[Unique["PolyLogLadder`O"], {Length[others]}];
  g = Expand[f /. Thread[sp[[All, 1]] -> lam] /. Thread[others -> osym]];
  If[! PolynomialQ[g, lam], Throw[f, PolyLogLadder`split]];
  terms = If[Head[g] === Plus, List @@ g, {g}];
  (* group by the monomial in the split kernels, partial fractions in x *)
  Do[
    With[{mono = Times @@ (lam^Exponent[t, lam])},
      If[mono === 1, out += t,
        Do[
          With[{p0 = PolyLogLadder`TermPole[a, x]},
            With[{pp = If[p0 === None || ! TrueQ[Im[p0] == 0], 0, p0]},
              out += a (mono /. Thread[lam -> (PolyLogLadder`SplitAt[#[[1, 1]], #[[2]], x, pp] & /@ sp)])]],
          {a, With[{ap = Apart[Together[t/mono], x]}, If[Head[ap] === Plus, List @@ ap, {ap}]]}]]],
    {t, terms}];
  out /. Thread[osym -> others]], PolyLogLadder`split];

(* {al, be, ga, de} of a Moebius M(x) = (al x + be)/(ga x + de), or $Failed. *)
PolyLogLadder`Moebius[m_, x_] := Module[{t = Together[m], nu, dn, c},
  nu = Numerator[t]; dn = Denominator[t];
  If[! PolynomialQ[nu, x] || ! PolynomialQ[dn, x] || Exponent[nu, x] > 1 || Exponent[dn, x] > 1,
    Return[$Failed]];
  c = {Coefficient[nu, x, 1], Coefficient[nu, x, 0], Coefficient[dn, x, 1], Coefficient[dn, x, 0]};
  If[PolyLogLadder`ZeroQ[c[[1]] c[[4]] - c[[2]] c[[3]]], $Failed, c]];

(* INT R(x) Lam^k dx, Lam a function with d/dx Lam = M'/M for the Moebius M
   (Lam = Log[M] itself, or a difference of two linear logs). *)
PolyLogLadder`KernelPiece[R_, k_, m_, lam_, x_] := Module[{c, u, xs, piece},
  c = PolyLogLadder`Moebius[m, x];
  If[c === $Failed, Return[$Failed]];
  u = Unique["PolyLogLadder`W"];
  xs = (c[[4]] u - c[[2]])/(c[[1]] - c[[3]] u);
  piece = PolyLogLadder`Ladder[Together[(R /. x -> xs) D[xs, u]], u, k];
  If[piece === $Failed, Return[$Failed]];
  (piece /. Log[u] -> lam) /. u -> m];

(* the same, oriented |u(37/100)| <= 1 when M is real (real coefficients):
   then u is real on the real line, and the ladder's Log[1 - u/c] and
   PolyLog[n, u/c] can only jump where u/c reaches 0 or 1 -- a singularity of
   the integrand itself.  A complex M (a split kernel 1 + I x) keeps u = M:
   1/(1 + I x) runs along an arc on which those pieces can cross a cut
   unseen. *)
PolyLogLadder`OrientedPiece[R_, k_, m_, lam_, x_] := Module[{mv, cs},
  cs = PolyLogLadder`Moebius[m, x];
  If[cs === $Failed || ! AllTrue[cs, NumericQ[#] && TrueQ[Im[#] == 0] &],
    Return[PolyLogLadder`KernelPiece[R, k, m, lam, x]]];
  mv = N[Abs[m /. x -> 37/100], 20];
  If[NumericQ[mv] && TrueQ[mv > 1],
    PolyLogLadder`KernelPiece[(-1)^k R, k, Together[1/m], -lam, x],
    PolyLogLadder`KernelPiece[R, k, m, lam, x]]];

(* Main entry: the log-kernel class. *)
PolyLogLadder`Run[f0_, x_Symbol] := Catch[Module[
  {f = f0, logs, Ls, g, rules, deg, res = 0, r0, piece, coeffs, chk, lin, cross, pr, wi, wj},
  logs = Union[Cases[f, Log[u_] /; ! FreeQ[u, x], {0, Infinity}]];
  If[logs === {}, Return[$Failed]];
  (* split the log of a real polynomial without real roots into linear logs *)
  If[AnyTrue[logs, PolynomialQ[#[[1]], x] && Exponent[#[[1]], x] >= 2 &],
    f = PolyLogLadder`SplitAll[f, x];
    logs = Union[Cases[f, Log[u_] /; ! FreeQ[u, x], {0, Infinity}]]];
  (* every log kernel Moebius in x (linear the common case) *)
  If[! AllTrue[logs, PolyLogLadder`Moebius[#[[1]], x] =!= $Failed &], Return[$Failed]];
  lin = (PolynomialQ[#[[1]], x] && Exponent[#[[1]], x] == 1) & /@ logs;
  Ls = Table[Unique["PolyLogLadder`L"], {Length[logs]}];
  rules = Thread[logs -> Ls];
  g = f /. rules;
  (* rational in x, polynomial in the kernels *)
  (* (constant logs such as the Log[P(p)] of a split are coefficients) *)
  If[! FreeQ[g /. Thread[Ls -> 1], (Log | Exp | Sin | Cos | Tan | Cot | Sec | Csc |
        Sinh | Cosh | Tanh | Coth | Sech | Csch | ArcTan | ArcSin | ArcCos | PolyLog)[___, v_, ___] /; ! FreeQ[v, x]] ||
     ! FreeQ[g, Power[E, v_] /; ! FreeQ[v, x]], Return[$Failed]];
  g = Expand[g];
  If[! AllTrue[Ls, PolynomialQ[g, #] &], Return[$Failed]];
  coeffs = {};
  r0 = g /. Thread[Ls -> 0];
  Do[
    deg = Exponent[g, Ls[[i]]];
    Do[
      AppendTo[coeffs, {i, kk, Coefficient[g /. Thread[Delete[Ls, i] -> 0], Ls[[i]], kk]}],
      {kk, 1, deg}],
    {i, Length[Ls]}];
  (* bilinear cross terms of two linear kernels; nothing of higher mixed degree *)
  cross = {};
  Do[
    If[lin[[i]] && lin[[j]],
      pr = Coefficient[Coefficient[g, Ls[[i]], 1], Ls[[j]], 1] /. Thread[Ls -> 0];
      If[! PolyLogLadder`ZeroQ[pr], AppendTo[cross, {i, j, pr}]]],
    {i, Length[Ls]}, {j, i + 1, Length[Ls]}];
  chk = Together[g - r0 - Total[(Ls[[#[[1]]]]^#[[2]] #[[3]]) & /@ coeffs] -
                 Total[(Ls[[#[[1]]]] Ls[[#[[2]]]] #[[3]]) & /@ cross]];
  If[! PolyLogLadder`ZeroQ[chk], Return[$Failed]];
  (* L_i L_j = (L_i^2 + L_j^2 - (L_i - L_j)^2)/2 *)
  Do[
    With[{i = c[[1]], j = c[[2]], R = c[[3]]},
      AppendTo[coeffs, {i, 2, R/2}];
      AppendTo[coeffs, {j, 2, R/2}];
      wi = logs[[i, 1]]; wj = logs[[j, 1]];
      piece = PolyLogLadder`OrientedPiece[-R/2, 2, Together[wi/wj], logs[[i]] - logs[[j]], x];
      If[piece === $Failed, Throw[$Failed, PolyLogLadder`fail]];
      res += piece],
    {c, cross}];
  Do[
    With[{i = c[[1]], kk = c[[2]], R = c[[3]]},
      piece = PolyLogLadder`OrientedPiece[R, kk, logs[[i, 1]], logs[[i]], x];
      If[piece === $Failed, Throw[$Failed, PolyLogLadder`fail]];
      res += piece],
    {c, coeffs}];
  If[! PolyLogLadder`ZeroQ[r0],
    piece = Integrate[r0, x];
    If[! FreeQ[piece, Integrate], Return[$Failed]];
    res += piece];
  res], PolyLogLadder`fail];

(* Trig route: polynomial x-weight times trig/hyperbolic of linear arguments.
   Every argument must be a rational multiple of one base frequency; the
   exponentials of TrigToExp are written as powers of the single kernel
   T = E^(base x), so the integrand is Sum_n x^n N_n(T)/Den(T).  Partial
   fractions in T over the roots of Den give three kinds of term:
     x^n T^k                    elementary (an exponential polynomial),
     x^n/(T - rho)              the exponential polylogarithm ladder
        INT = Sum_k -(1/rho) n!/(n-k)! / base^(k+1) x^(n-k) PolyLog[k+1, rho/T],
     x^n/(T - rho)^m, m >= 2    (1/(m-1)!) d^(m-1)/drho^(m-1) of the m = 1 form,
   each an exact identity (re-checked by the C stage by differentiation). *)
PolyLogLadder`TrigHeads = Sin | Cos | Tan | Cot | Sec | Csc | Sinh | Cosh | Tanh | Coth | Sech | Csch;

(* INT x^n/(E^(base x) - rho) dx, rho != 0 (symbolic rho allowed). *)
PolyLogLadder`ExpLadder[n_, rho_, base_, x_] :=
  Sum[-(1/rho) n!/(n - k)!/base^(k + 1) x^(n - k) PolyLog[k + 1, rho/E^(base x)], {k, 0, n}];

PolyLogLadder`TrigRun[f_, x_Symbol] := Module[
  {args, cs, c0, ratios, q, base, g, T, num, den, res = 0, nmax, Nn, poly, rem, roots, qq, h, cb, rs, F1, term, sub},
  If[FreeQ[f, PolyLogLadder`TrigHeads], Return[$Failed]];
  If[! FreeQ[f, Log | PolyLog], Return[$Failed]];
  args = Union[Cases[f, PolyLogLadder`TrigHeads[u_] :> u, {0, Infinity}]];
  If[! AllTrue[args, PolynomialQ[#, x] && Exponent[#, x] == 1 && Coefficient[#, x, 0] === 0 &], Return[$Failed]];
  (* an explicit polynomial weight in x: a pure trig rational belongs to Weierstrass *)
  If[FreeQ[f /. PolyLogLadder`TrigHeads[_] -> 1, x], Return[$Failed]];
  If[! FreeQ[f, Sin | Cos | Tan | Cot | Sec | Csc] && ! FreeQ[f, Sinh | Cosh | Tanh | Coth | Sech | Csch], Return[$Failed]];
  cs = Coefficient[#, x, 1] & /@ args;
  c0 = First[cs];
  ratios = Simplify[cs/c0];
  If[! AllTrue[ratios, Head[#] === Integer || Head[#] === Rational &], Return[$Failed]];
  q = LCM @@ (Denominator /@ ratios);
  base = If[FreeQ[f, Sin | Cos | Tan | Cot | Sec | Csc], c0/q, I c0/q];
  T = Unique["PolyLogLadder`T"];
  g = TrigToExp[f] /. Power[E, z_] /; ! FreeQ[z, x] :> T^Simplify[z/(base x)];
  If[! FreeQ[Cases[g, Power[E, z_] :> z, {0, Infinity}], x], Return[$Failed]];
  g = Together[g];
  num = Expand[Numerator[g]]; den = Expand[Denominator[g]];
  If[! FreeQ[den, x] || ! PolynomialQ[num, x], Return[$Failed]];
  If[! PolynomialQ[num, T] || ! PolynomialQ[den, T], Return[$Failed]];
  (* pull the T^j factor of the denominator out as negative powers *)
  nmax = Exponent[num, x];
  roots = Select[PolyLogLadder`RootsMult[den, T] /. $Failed -> {}, ! PolyLogLadder`ZeroQ[#[[1]]] &];
  If[PolyLogLadder`RootsMult[den, T] === $Failed, Return[$Failed]];
  rs = Unique["PolyLogLadder`rho"];
  Do[
    Nn = Coefficient[num, x, n];
    If[Nn === 0, Continue[]];
    (* exponential-polynomial part: Laurent terms in T around 0 and the
       polynomial quotient, via the partial-fraction remainder *)
    poly = PolynomialQuotient[Nn, den, T];
    rem = PolynomialRemainder[Nn, den, T];
    res += Integrate[Expand[x^n (poly /. T -> E^(base x))], x];
    Do[
      With[{b = rt[[1]], m = rt[[2]]},
        qq = PolynomialQuotient[den, Expand[(T - b)^m], T];
        h = rem/qq;
        Do[
          cb = Simplify[(D[h, {T, m - j}] /. T -> b)/(m - j)!];
          If[! PolyLogLadder`ZeroQ[cb],
            F1 = PolyLogLadder`ExpLadder[n, rs, base, x];
            term = If[j == 1, F1, D[F1, {rs, j - 1}]/(j - 1)!] /. rs -> b;
            res += cb term],
          {j, 1, m}]],
      {rt, roots}];
    (* a root T = 0 of the denominator: the remainder over T^j T-part *)
    With[{z = Select[PolyLogLadder`RootsMult[den, T], PolyLogLadder`ZeroQ[#[[1]]] &]},
      If[z =!= {},
        With[{m0 = z[[1, 2]]},
          qq = PolynomialQuotient[den, T^m0, T];
          h = rem/qq;
          Do[
            cb = Simplify[(D[h, {T, m0 - j}] /. T -> 0)/(m0 - j)!];
            If[! PolyLogLadder`ZeroQ[cb],
              res += cb Integrate[x^n E^(-j base x), x]],
            {j, 1, m0}]]]],
    {n, 0, nmax}];
  If[! FreeQ[res, Integrate], Return[$Failed]];
  res];


(* Log of a trig rational: P(x) Log[G], P a polynomial, G a rational function of
   trig/hyperbolic functions of linear arguments with one base frequency.  With
   T = E^(base x) and G = K T^e Prod (T - rho)^m, d/dx Log[G] equals
   d/dx [e base x + Sum m Log[1 - T/rho]] exactly, so
     Log[G] = D(x) + L(x),   L = e base x + Sum m Log[1 - T/rho],
   with D = Log[G] - L piecewise constant (it jumps only where a factor crosses
   a branch cut).  Hence
     INT P Log[G] = D(x) INT P + INT P L,
   exact wherever D is continuous -- D is kept symbolic, never sampled -- and
   INT x^n Log[1 - T/rho] = -INT x^n PolyLog[1, T/rho] is the exponential ladder
     INT x^n PolyLog[s, T/rho] = Sum (-1)^j n!/(n-j)! x^(n-j) PolyLog[s+1+j, T/rho]/base^(j+1).
   A definite integral across a jump of D is caught by Newton-Leibniz's
   continuity cross-check (it fires on PolyLog/Log of an exponential). *)
(* INT x^n PolyLog[s, T/rho] dx, T = E^(base x) *)
PolyLogLadder`LTRoots[p_, T_] := If[FreeQ[p, T], {}, PolyLogLadder`RootsMult[p, T]];
PolyLogLadder`LTLiLadder[n_, s_, rho_, base_, x_] :=
  Sum[(-1)^j n!/(n - j)! x^(n - j) PolyLog[s + 1 + j, E^(base x)/rho]/base^(j + 1), {j, 0, n}];
PolyLogLadder`LTRun[f_, x_Symbol] := Catch[Module[{logs, Ls, g, P, G, args, cs, c0, ratios, q, base, T, GT, num, den, e, roots, L, LI, res = 0, Q, dexact},
  logs = Union[Cases[f, Log[u_] /; ! FreeQ[u, x], {0, Infinity}]];
  If[Length[logs] != 1, Throw[$Failed, PolyLogLadder`fail]];
  G = logs[[1, 1]];
  If[FreeQ[G, PolyLogLadder`TrigHeads], Throw[$Failed, PolyLogLadder`fail]];
  P = Together[f/logs[[1]]];
  If[! PolynomialQ[P, x], Throw[$Failed, PolyLogLadder`fail]];
  args = Union[Cases[G, PolyLogLadder`TrigHeads[u_] :> u, {0, Infinity}]];
  If[! AllTrue[args, PolynomialQ[#, x] && Exponent[#, x] == 1 && Coefficient[#, x, 0] === 0 &], Throw[$Failed, PolyLogLadder`fail]];
  If[! FreeQ[G /. PolyLogLadder`TrigHeads[_] -> 1, x], Throw[$Failed, PolyLogLadder`fail]];
  cs = Coefficient[#, x, 1] & /@ args; c0 = First[cs];
  ratios = Simplify[cs/c0];
  If[! AllTrue[ratios, Head[#] === Integer || Head[#] === Rational &], Throw[$Failed, PolyLogLadder`fail]];
  q = LCM @@ (Denominator /@ ratios);
  base = If[FreeQ[G, Sin | Cos | Tan | Cot | Sec | Csc], c0/q, I c0/q];
  T = Unique["PolyLogLadder`T"];
  GT = Together[TrigToExp[G] /. Power[E, z_] /; ! FreeQ[z, x] :> T^Simplify[z/(base x)]];
  If[! FreeQ[GT, x], Throw[$Failed, PolyLogLadder`fail]];
  num = Numerator[GT]; den = Denominator[GT];
  (* e: net power of T as T -> 0 *)
  e = Exponent[num, T, Min] - Exponent[den, T, Min];
  roots = Join[{#[[1]], #[[2]]} & /@ PolyLogLadder`LTRoots[Expand[num/T^Exponent[num, T, Min]], T],
               {#[[1]], -#[[2]]} & /@ PolyLogLadder`LTRoots[Expand[den/T^Exponent[den, T, Min]], T]];
  If[! FreeQ[roots, $Failed], Throw[$Failed, PolyLogLadder`fail]];
  L = e base x + Total[(#[[2]] Log[1 - E^(base x)/#[[1]]]) & /@ roots];
  Q = Integrate[P, x];
  (* INT P L dx *)
  res = Integrate[Expand[P e base x], x];
  Do[
    With[{c = Coefficient[Expand[P], x, n]},
      If[c =!= 0,
        res += c Total[(-#[[2]] PolyLogLadder`LTLiLadder[n, 1, #[[1]], base, x]) & /@ roots]]],
    {n, 0, Exponent[P, x]}];
  (* D = Log[G] - L is piecewise constant: G = K T^e Prod (T - rho)^m and
     (T - rho) = -rho (1 - T/rho) give D = Log[K Prod (-rho)^m] + 2 Pi I n.  The
     integer n is read numerically at x = 37/100 (where the acceptance check
     also differentiates); a symbolic G keeps D as Log[G] - L. *)
  dexact = Module[{K, C0, Dn, n},
    K = Coefficient[num, T, Exponent[num, T]]/Coefficient[den, T, Exponent[den, T]];
    C0 = K Times @@ ((-#[[1]])^#[[2]] & /@ roots);
    Dn = Quiet[N[((Log[G] - L) /. x -> 37/100) - Log[C0], 20]];   (* Re[Dn] is a zero *)
    If[! NumericQ[Dn] || Abs[Re[Dn]] > 10^-8, Log[G] - L,
      n = Round[Im[Dn]/(2 Pi)];
      If[Abs[Im[Dn] - 2 Pi n] > 10^-8, Log[G] - L, Log[C0] + 2 Pi I n]]];
  dexact Q + res], PolyLogLadder`fail];


(* Acceptance check: D[F, x] - f vanishes at two points, every other free
   symbol fixed at a generic rational.  A differentiation check, not quadrature;
   the constructions above are exact, so this only guards against an
   implementation slip. *)
PolyLogLadder`Check[F_, f_, x_Symbol] := Module[{syms, inst, d, ok = True},
  (* parameters: the non-numeric atoms other than x (Pi, E, Catalan, ... are
     NumericQ; Context[] is not usable to separate system symbols here) *)
  syms = DeleteCases[Union[Cases[{F, f}, s_Symbol /; ! NumericQ[s] && s =!= x, {0, Infinity}]], x];
  inst = Thread[syms -> Take[{23/29, 31/37, 13/17, 41/43, 7/11, 19/23, 29/31, 3/7}, UpTo[Length[syms]]]];
  If[Length[inst] < Length[syms], Return[False]];
  Do[
    d = Quiet[N[(D[F, x] - f) /. inst /. x -> pt, 30]];   (* a zero residual: N::meprec expected *)
    If[! NumericQ[d] || Abs[d] > 10^-12 (1 + Abs[N[f /. inst /. x -> pt, 30]]), ok = False],
    {pt, {37/100, 71/100}}];
  ok];

PolyLogLadder`Integrate[f_, x_Symbol] := Module[{r},
  r = PolyLogLadder`Run[f, x];
  If[r === $Failed, r = PolyLogLadder`PLRun[f, x]];
  If[r === $Failed, r = PolyLogLadder`TrigRun[f, x]];
  If[r === $Failed, r = PolyLogLadder`LTRun[f, x]];
  If[r === $Failed || ! FreeQ[r, Integrate], Return[$Failed]];
  If[PolyLogLadder`Check[r, f, x], r, $Failed]];


(* ---------------------------------------------------------------------------
   POLYLOG ROUTE.  f = R0(x) + Sum_i R_i(x) PolyLog[n_i, a_i(x)], every R
   rational in x, every n_i an integer >= 2 and every a_i either

     Moebius      a = (al x + be)/(ga x + de), al de - be ga != 0, or
     monomial     a = c x^k, k an integer, |k| >= 2.

   METHOD.  Integration by parts lowers the polylogarithm's order one step at a
   time, d/dt Li_n(w) = Li_(n-1)(w) w'/w, down to Li_1(w) = -Log[1 - w], where
   the log-kernel ladder above (or the cascade) finishes.  To keep w'/w simple
   the Moebius case first substitutes u = a(x) (x = (de u - be)/(al - ga u) is
   Moebius again, so the integrand stays rational: Q(u) = R(x(u)) x'(u)); then
   w = u and w'/w = 1/u.  A monomial argument keeps t = x, with w'/w = k/x.
   In both cases, with Q = r/t + s/(t - 1) + Q~ (partial fractions in t) and
   P = INT Q~ dt:

     INT Q Li_n(w) dt = (r/k) Li_(n+1)(w) + P Li_n(w) - k INT (P/t) Li_(n-1)(w) dt
                        + s U_n(t),

   which needs P rational: every pole of Q other than t = 0 (and t = 1 when
   w = t) must have zero residue.  It always does when R is a polynomial and
   the argument Moebius (x^m dx is exact); a nonzero residue at another point
   leaves INT Li_n(w)/(t - c) dt, which is not a classical polylogarithm, and
   the route declines.  The t = 1 pole (w = t only) closes for n <= 3:

     U_n = INT Li_n(t)/(t - 1) dt
     U_2 = Log[1 - t] Li_2(t) + INT Log[1 - t]^2/t dt        (the log ladder)
     U_3 = Log[1 - t] Li_3(t) + Li_2(t)^2/2.

   Each identity differentiates back exactly; the stage re-checks the result
   by differentiation regardless.
   --------------------------------------------------------------------------- *)

PolyLogLadder`plSavedContextPath = $ContextPath;
$ContextPath = {"System`"};
Begin["PolyLogLadder`Private`"];

PolyLogLadder`PLRatQ[e_, t_] := With[{tg = Together[e]},
  PolynomialQ[Numerator[tg], t] && PolynomialQ[Denominator[tg], t]];

(* U_n(t) = INT PolyLog[n, t]/(t - 1) dt, or $Failed. *)
PolyLogLadder`PLUnitPole[2, t_] := With[{sub = PolyLogLadder`Run[Log[1 - t]^2/t, t]},
  If[sub === $Failed, $Failed, Log[1 - t] PolyLog[2, t] + sub]];
PolyLogLadder`PLUnitPole[3, t_] := Log[1 - t] PolyLog[3, t] + PolyLog[2, t]^2/2;
PolyLogLadder`PLUnitPole[_, _] := $Failed;

(* INT Q(t) PolyLog[n, w] dt, w = t (k = 1) or w = c t^k: t w'/w = k. *)
PolyLogLadder`PLCore[Q_, n_Integer, w_, t_Symbol, k_Integer] := Catch[Module[
  {qq = Together[Q], ap, terms, r = 0, s = 0, pp = 0, pt, sub, us, res},
  If[qq === 0, Return[0]];
  If[n == 1,
    res = If[w === t, PolyLogLadder`Run[-qq Log[1 - t], t], Integrate[-qq Log[1 - w], t]];
    If[res === $Failed || ! FreeQ[res, Integrate], Throw[$Failed, PolyLogLadder`plfail]];
    Return[res]];
  ap = Apart[qq, t];
  terms = If[Head[ap] === Plus, List @@ ap, {ap}];
  Do[
    Which[
      FreeQ[Together[tm t], t], r += Together[tm t],
      w === t && FreeQ[Together[tm (t - 1)], t], s += Together[tm (t - 1)],
      True,
        pt = Integrate[tm, t];
        If[! FreeQ[pt, Integrate] || ! PolyLogLadder`PLRatQ[pt, t], Throw[$Failed, PolyLogLadder`plfail]];
        pp += pt],
    {tm, terms}];
  pp = Together[pp];
  sub = PolyLogLadder`PLCore[pp/t, n - 1, w, t, k];
  If[sub === $Failed, Throw[$Failed, PolyLogLadder`plfail]];
  res = r/k PolyLog[n + 1, w] + pp PolyLog[n, w] - k sub;
  If[s =!= 0,
    us = PolyLogLadder`PLUnitPole[n, t];
    If[us === $Failed, Throw[$Failed, PolyLogLadder`plfail]];
    res += s us];
  res], PolyLogLadder`plfail];

(* INT R(x) PolyLog[n, a] dx for one kernel, or $Failed. *)
PolyLogLadder`PLKernel[R_, n_Integer, a_, x_Symbol] := Module[
  {at = Together[a], nu, dn, al, be, ga, d0, u, xs, qu, g, kk},
  nu = Numerator[at]; dn = Denominator[at];
  If[PolynomialQ[nu, x] && PolynomialQ[dn, x] && Exponent[nu, x] <= 1 && Exponent[dn, x] <= 1,
    al = Coefficient[nu, x, 1]; be = Coefficient[nu, x, 0];
    ga = Coefficient[dn, x, 1]; d0 = Coefficient[dn, x, 0];
    If[PolyLogLadder`ZeroQ[al d0 - be ga], Return[$Failed]];
    u = Unique["PolyLogLadder`u"];
    xs = (d0 u - be)/(al - ga u);
    qu = Together[(R /. x -> xs) D[xs, u]];
    g = PolyLogLadder`PLCore[qu, n, u, u, 1];
    Return[If[g === $Failed, $Failed, g /. u -> a]]];
  (* monomial c x^k *)
  kk = Together[x D[a, x]/a];
  If[! IntegerQ[kk] || Abs[kk] < 2 || ! FreeQ[Together[a/x^kk], x], Return[$Failed]];
  PolyLogLadder`PLCore[R, n, a, x, kk]];

PolyLogLadder`PLRun[f_, x_Symbol] := Catch[Module[{pls, ks, g, r0, res = 0, ri, piece},
  pls = Union[Cases[f, PolyLog[_, v_] /; ! FreeQ[v, x], {0, Infinity}]];
  If[pls === {}, Throw[$Failed, PolyLogLadder`plfail]];
  If[! AllTrue[pls, IntegerQ[#[[1]]] && #[[1]] >= 2 &], Throw[$Failed, PolyLogLadder`plfail]];
  ks = Table[Unique["PolyLogLadder`K"], {Length[pls]}];
  g = f /. Thread[pls -> ks];
  If[! FreeQ[g, PolyLog], Throw[$Failed, PolyLogLadder`plfail]];
  r0 = g /. Thread[ks -> 0];
  Do[
    ri = D[g, ks[[i]]];
    If[! FreeQ[ri, Alternatives @@ ks] || ! PolyLogLadder`PLRatQ[ri, x], Throw[$Failed, PolyLogLadder`plfail]];
    piece = PolyLogLadder`PLKernel[ri, pls[[i, 1]], pls[[i, 2]], x];
    If[piece === $Failed, Throw[$Failed, PolyLogLadder`plfail]];
    res += piece,
    {i, Length[ks]}];
  (* linear in the kernels, rational remainder *)
  If[! PolyLogLadder`ZeroQ[g - r0 - Total[Table[D[g, ks[[i]]] ks[[i]], {i, Length[ks]}]]],
    Throw[$Failed, PolyLogLadder`plfail]];
  If[! PolyLogLadder`ZeroQ[r0],
    If[! PolyLogLadder`PLRatQ[r0, x], Throw[$Failed, PolyLogLadder`plfail]];
    piece = Integrate[r0, x];
    If[! FreeQ[piece, Integrate], Throw[$Failed, PolyLogLadder`plfail]];
    res += piece];
  PolyLogLadder`PLTidy[res, x]], PolyLogLadder`plfail];

(* Presentation only, and every step an identity on the whole complex plane:
   ArcTanh[v] = (Log[1 + v] - Log[1 - v])/2 (the principal branches agree);
   kernel arguments over one denominator (the same number, so no branch moves);
   a positive rational content pulled out of a Log, Log[c u] = Log[c] + Log[u]
   for c > 0; then the antiderivative collected in its kernels with factored
   coefficients.  The Log form lets the singular pieces at an endpoint cancel
   term by term (c(x) Log[1 - x] with c(1) = 0), where ArcTanh[1 - 2 x] next to
   Log[1 - x] leaves the cancellation to a limit of two divergent kernels. *)
PolyLogLadder`PLLogSplit[v_, x_] := Module[{tv = Together[v], nu, dn, cn, cd, k},
  nu = Numerator[tv]; dn = Denominator[tv];
  If[! PolynomialQ[nu, x] || ! PolynomialQ[dn, x], Return[Log[tv]]];
  cn = CoefficientList[nu, x]; cd = CoefficientList[dn, x];
  If[! AllTrue[Join[cn, cd], IntegerQ], Return[Log[tv]]];
  k = (GCD @@ cn)/(GCD @@ cd);
  (* keep the sign with the polynomial: only a positive content splits off *)
  If[k === 1 || ! TrueQ[k > 0], Log[tv], Log[k] + Log[Together[tv/k]]]];

PolyLogLadder`PLTidy[e_, x_Symbol] := Module[{g, ks},
  g = e /. ArcTanh[v_] /; ! FreeQ[v, x] :> (Log[1 + v] - Log[1 - v])/2;
  g = g /. {Log[v_] /; ! FreeQ[v, x] :> PolyLogLadder`PLLogSplit[v, x],
            PolyLog[m_, v_] :> PolyLog[m, Together[v]]};
  ks = Union[Cases[g, kn : (Log | PolyLog | ArcTan)[__] /; ! FreeQ[kn, x], {0, Infinity}]];
  If[ks === {}, Together[g], Collect[Expand[g], ks, Factor]]];

End[];
$ContextPath = PolyLogLadder`plSavedContextPath;
