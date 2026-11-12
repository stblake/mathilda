(* families.m -- integrand families shared by the Bronstein stress harness
   (tools/bronstein_stress/gen.m) and the integration shootout
   (tools/integration_shootout/gen.m).  Every family is a Function[{th, n}, F]
   returning an ELEMENTARY F built from the monomial th and a radical degree n;
   the integrand is D[F, x].  The helpers draw from the global random state,
   which the caller seeds; nothing here reseeds.

   mixedFamilies: Shape A (a radical whose radicand involves the top Log/Exp
   monomial).  theta[] draws a random monomial, lin[th] a radicand linear in it.
   The order of mixedFamilies is load-bearing: the Bronstein nightly's seeds and
   its --max-open baseline were recorded with it. *)
r[a_, b_] := RandomInteger[{a, b}];
nz[a_, b_] := RandomChoice[DeleteCases[Range[a, b], 0]];
theta[] := RandomChoice[{Log[x], Log[x + r[1, 3]], Log[nz[1, 3] x + r[1, 3]], E^x, E^(nz[-2, 2] x),
   E^(x^2), Log[1 + E^x], Log[x^2 + 1]}];
lin[th_] := nz[1, 3] th + r[-2, 2] x + nz[-3, 3];
mixedFamilies = {
  (* algebraic part only *)
  "alg" -> Function[{th, n}, With[{Y = lin[th]^(1/n)}, (r[1, 3] x + r[-2, 2] th + r[-1, 1]) Y^RandomChoice[{1, -1, n - 1}]]],
  (* log of x + Y: residues at finite places *)
  "logxy" -> Function[{th, n}, With[{Y = lin[th]^(1/n)}, r[0, 2] x Y + nz[-2, 2] Log[x + r[1, 2] + Y]]],
  (* arcsinh: residues at the special places *)
  "asinh" -> Function[{th, n}, With[{Y = Sqrt[th^2 + nz[1, 3]]}, nz[-2, 2] Log[th + Y] + r[0, 1] Y]],
  (* genus 1: Log[p + Y] on Y^2 = th^3 + c *)
  "genus1" -> Function[{th, n}, With[{Y = Sqrt[th^3 + nz[1, 3]]}, Log[r[1, 2] + Y] + r[0, 1] th Y]],
  (* algebraic residues: Log of a quadratic in Y over Q(sqrt) *)
  "algres" -> Function[{th, n}, With[{Y = Sqrt[th^2 + nz[1, 2]]}, Log[(th + 1 + Y)/(th - 1 + Y)] + r[0, 1] Y/th]],
  (* rational multiple of a component: b(th) Y with b rational in th *)
  "ratb" -> Function[{th, n}, With[{Y = lin[th]^(1/n)}, Y/(th + nz[1, 3]) + r[0, 1] Y^(n - 1)]],
  (* BRONSTEIN_PLAN.md gap #3: the shapes the six families above never produce *)
  (* constant factors inside the kernels (Log[2 x + 2], E^(x + 1), 3 E^x) *)
  "constk" -> Function[{th, n}, With[{t = RandomChoice[{Log[nz[2, 3] x], Log[nz[2, 3] x + nz[2, 3]], E^(x + nz[1, 2]), nz[2, 3] E^x, Log[nz[2, 3] (x^2 + 1)]}]},
    With[{Y = lin[t]^(1/n)}, r[1, 2] x Y + nz[-2, 2] Log[x + r[1, 2] + Y]]]],
  (* the radicand on two monomials of a two-level tower (both stay in k(theta)) *)
  "tower2" -> Function[{th, n}, With[{Y = (nz[1, 2] E^x + nz[1, 2] Log[x] + r[-1, 2])^(1/n)},
    RandomChoice[{r[1, 2] x Y, Log[x + Y], Y/(Log[x] + nz[1, 3]), Log[1 + Y] + Y}]]],
  (* non-constant leading coefficient and content in k: x theta + c, x (theta + c) *)
  "lccontent" -> Function[{th, n}, With[{Y = RandomChoice[{(x th + nz[1, 3])^(1/n), (x (th + nz[1, 3]))^(1/n), ((x + 1) th^2 + 1)^(1/n)}]},
    RandomChoice[{th Y, Y/x, x Y + Log[x + Y], Log[1 + Y]}]]],
  (* exponential, cube root: logarithms with residues at the special places *)
  "expcube" -> Function[{th, n}, With[{Y = (RandomChoice[{E^x, E^(2 x), E^(x^2)}] + nz[1, 3])^(1/3)},
    nz[-2, 2] Log[Y + RandomChoice[{-1, 1, 2}]] + r[0, 1] Y]]
};
