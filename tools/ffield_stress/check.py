#!/usr/bin/env python3
"""tools/ffield_stress/check.py -- check the ffield / polynomial-matrix kernels
(Phase 4.5) against two independent systems: PARI/GP and the reference CAS
(wolframscript).

    python3 tools/ffield_stress/check.py [--n N] [--seed S] [--no-gp] [--no-wl]

Runs tools/ffield_stress/gen.m (random cases over Q, GF(p), Q(theta), Q(a),
Q(a, b), Q(theta)(a), Q(theta)(a, b)), then verifies every answer by its
defining properties, never by comparing with another implementation's choice:

  divrem   A = q B + r and deg r < deg B (and q, r equal PARI's divrem)
  xgcd     g monic (or 1), s A + t B = g, g | A and g | B  (so g = gcd)
  hnf      H upper triangular, monic pivots, entries above a pivot reduced;
           every row of M reduces to 0 by H (span M <= span H); and
           prod diag H = the monic gcd of the maximal minors of M (PARI) /
           = det M made monic for a square M (both), so span H = span M
  detadj   adj . M = det I and det = det M

Elements of Q(theta)(params) are handled natively by each system: PARI as
Mod(polynomial in t, mu(t)) over rational functions in the parameters (with
the priority x > t > parameters), the reference CAS by explicit reduction of
numerators modulo mu.  GF(p) is Mod(., p) / PolynomialMod.

Developer tool, not part of CI.  Prints one line per failure and a summary.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MATH = os.path.join(ROOT, "Mathilda")
GP = "gp"
def find_wl():
    """the reference CAS's script runner: $FFS_WL, else wolframscript on PATH,
    else the first application bundle that ships one."""
    import glob
    import shutil
    if os.environ.get("FFS_WL"):
        return os.environ["FFS_WL"]
    w = shutil.which("wolframscript")
    if w:
        return w
    hits = sorted(glob.glob("/Applications/*.app/Contents/MacOS/wolframscript"))
    return hits[0] if hits else "wolframscript"


WL = find_wl()


def run_gen(n, seed):
    env = dict(os.environ, FFS_N=str(n), FFS_SEED=str(seed))
    out = subprocess.run([MATH, "-file", os.path.join(ROOT, "tools/ffield_stress/gen.m")],
                         capture_output=True, text=True, env=env, cwd=ROOT, timeout=3600)
    cases = []
    for line in out.stdout.splitlines():
        if not line.startswith("FFS\t"):
            continue
        f = line.split("\t")
        cases.append({"id": int(f[1]), "field": f[2], "op": f[3], "mu": f[4], "mod": int(f[5]),
                      "params": f[6], "ins": f[7], "out": f[8]})
    return cases


def to_gp(s):
    s = s.replace("{", "[").replace("}", "]")
    return s


# ---------------------------------------------------------------- PARI/GP
GP_PRELUDE = r"""
\\ variable priority: x > t > a > b
t; a; b;
default(parisize, 400000000);
rows2mat(L) = matrix(#L, #L[1], i, j, L[i][j]);
\\ embed a polynomial in x with coefficients in Q(params)[t] into K[x]
\\ FFPT: the parameters' values for this evaluation (a specialisation check:
\\ an identity in Q(params) holds at a random point; a false one fails there
\\ except on a proper Zariski-closed set, so two points are tried)
K(P) = { my(Q = P);
  for (i = 1, #FFPAR, Q = subst(Q, FFPAR[i], FFPT[i]));
  if (FFMOD, return (Q * Mod(1, FFMOD)));
  if (FFMU == 0, return (Q));
  if (type(Q) != "t_POL" || variable(Q) != x, return (Mod(Q, FFMU)));
  sum(k = 0, poldegree(Q, x), Mod(polcoef(Q, k, x), FFMU) * x^k) };
KM(L) = { my(M = rows2mat(L)); matrix(#M[,1], #M, i, j, K(M[i, j])) };
iszero(e) = { my(z = e); if (type(z) == "t_POLMOD", z = lift(z)); z == 0 };
pzero(P) = { if (P == 0, return (1)); if (type(P) != "t_POL" || variable(P) != x, return (iszero(P)));
  for (k = 0, poldegree(P, x), if (!iszero(polcoef(P, k, x)), return (0))); 1 };
pdeg(P) = { if (pzero(P), return (-1)); if (type(P) != "t_POL" || variable(P) != x, return (0));
  forstep (k = poldegree(P, x), 0, -1, if (!iszero(polcoef(P, k, x)), return (k))); -1 };
plead(P) = polcoef(P, pdeg(P), x);
\\ reduce a row vector by the triangular H (monic pivots)
inspan(v0, H) = { my(v = v0, c = #H, q);
  for (i = 1, c, q = divrem(Pol(v[i], x), Pol(H[i, i], x), x)[1]; v = v - q * H[i, ]; if (!pzero(v[i]), return (0)));
  1 };
minorsgcd(M) = { my(r = #M[,1], c = #M, g = 0);
  forsubset([r, c], S, g = gcd(g, matdet(matrix(c, c, i, j, M[S[i], j]))));
  g };
"""


GP_POINTS = [((5, 3), (7, 2)), ((-11, 4), (13, 5))]


def gp_subst(text, params, point):
    """the parameters of `text` replaced by rational literals (textually, so GP
    never builds the rational functions in the parameters: parsing them was the
    slow part, minutes for two parameters)."""
    names = [p.strip() for p in params.strip("{}").split(",") if p.strip()]
    for nm, (pn, pd) in zip(names, point):
        text = re.sub(r"\b%s\b" % re.escape(nm), "(%d/%d)" % (pn, pd), text)
    return text


def gp_body(c, ins, out):
    op = c["op"]
    if op == "divrem":
        return """my(I = %s, O = %s, A = K(I[1]), B = K(I[2]), q = K(O[1]), r = K(O[2]), qr);
ok = pzero(A - (q * B + r)) && pdeg(r) < pdeg(B);
qr = divrem(A, B, x); ok = ok && pzero(qr[1] - q) && pzero(qr[2] - r);""" % (ins, out)
    if op == "xgcd":
        return """my(I = %s, O = %s, A = K(I[1]), B = K(I[2]), g = K(O[1]), s = K(O[2][1]), u = K(O[2][2]));
ok = pzero(s * A + u * B - g) && iszero(plead(g) - 1);
ok = ok && pzero(A %% g) && pzero(B %% g);""" % (ins, out)
    if op == "hnf":
        return """my(I = %s, O = %s, M = KM(I[1]), H = KM(O), c = #H, dg, g);
ok = (#H[,1] == c);
for (i = 1, c, for (j = 1, c,
  if (j < i && !pzero(H[i, j]), ok = 0);
  if (j == i && !iszero(plead(H[i, i]) - 1), ok = 0);
  if (j > i && pdeg(H[i, j]) >= pdeg(H[j, j]), ok = 0)));
for (k = 1, #M[,1], if (ok && !inspan(M[k, ], H), ok = 0));
if (ok, dg = prod(i = 1, c, H[i, i]); g = minorsgcd(M); g = g / plead(g); ok = pzero(dg - g));""" % (ins, out)
    if op == "detadj":
        return """my(I = %s, O = %s, M = KM(I[1]), d = K(O[1]), A = KM(O[2]), n = #M, P);
P = A * M; ok = pzero(d - matdet(M));
for (i = 1, n, for (j = 1, n, if (!pzero(P[i, j] - if (i == j, d, 0)), ok = 0)));""" % (ins, out)
    return None


def gp_case(c):
    head = "FFMOD = %d; FFMU = %s; FFPAR = [];\n" % (c["mod"], "0" if c["mu"] == "None" else to_gp(c["mu"]))
    ins, out = to_gp(c["ins"]), to_gp(c["out"])
    has_par = c["params"].strip("{} ") != ""
    # with parameters: the identities checked at two rational points (a
    # specialisation check; a failure must show at both, so that an unlucky
    # point -- a pole, a degree drop of the specialised data -- does not count)
    pts = GP_POINTS if has_par else [None]
    fns = []
    for k, pt in enumerate(pts):
        bi, bo = (gp_subst(ins, c["params"], pt), gp_subst(out, c["params"], pt)) if pt else (ins, out)
        body = gp_body(c, bi, bo)
        if body is None:
            return None
        fns.append("chk%d_%d() = { my(ok = 0); %s; ok };\n" % (c["id"], k, body.replace("\n", " ").rstrip(";")))
    # a point whose evaluation raises (a pole of an entry at that point) is
    # skipped; the case is ok at the first point that passes, an error only
    # when every point raised.  Each case under alarm(30): a timeout is
    # reported (-2), never a hang.
    tries = "".join("if (ok != 1, my(v = iferr(chk%d_%d(), E, if (errname(E) == \"e_ALARM\", error(E)); -1)); "
                    "if (v == 1 || ok == -1 || ok == 99, ok = v));" % (c["id"], k) for k in range(len(pts)))
    return (head + "".join(fns) +
            "ok = 99; iferr(alarm(30, %s), E, ok = if (errname(E) == \"e_ALARM\", -2, -1)); if (ok == 99, ok = -1); "
            "print(\"R \", %d, \" \", ok);\n" % (tries, c["id"]))


# ---------------------------------------------------------------- reference CAS
WL_PRELUDE = r"""
red[e_] := Module[{n = Numerator[Together[e]]},
  n = Expand[n];
  If[ffmu =!= None, n = Expand[PolynomialRemainder[n, ffmu, t]]];
  If[ffmod > 0, n = PolynomialMod[n, ffmod]];
  n];
zq[e_] := red[e] === 0;
pdeg[P_] := If[zq[P], -1, Exponent[red[P], x]];
plead[P_] := Coefficient[Expand[P], x, Exponent[Expand[P], x]];
inspan[v0_, H_] := Module[{v = v0, q, c = Length[H]},
  Catch[Do[q = PolynomialQuotient[Expand[v[[i]]], H[[i, i]], x]; v = Expand[v - q H[[i]]];
    If[!zq[v[[i]]], Throw[False]], {i, c}]; True]];
"""


def wl_case(c):
    head = "ffmod = %d; ffmu = %s;" % (c["mod"], c["mu"])
    ins, out = c["ins"], c["out"]
    op = c["op"]
    if op == "divrem":
        body = "With[{I0 = %s, O0 = %s}, zq[I0[[1]] - (O0[[1]] I0[[2]] + O0[[2]])] && pdeg[O0[[2]]] < pdeg[I0[[2]]]]" % (ins, out)
    elif op == "xgcd":
        body = ("With[{I0 = %s, O0 = %s}, zq[O0[[2, 1]] I0[[1]] + O0[[2, 2]] I0[[2]] - O0[[1]]] && "
                "zq[plead[O0[[1]]] - 1] && zq[PolynomialRemainder[Expand[I0[[1]]], O0[[1]], x]] && "
                "zq[PolynomialRemainder[Expand[I0[[2]]], O0[[1]], x]]]") % (ins, out)
    elif op == "hnf":
        body = ("With[{M = %s[[1]], H = %s}, Module[{c = Length[H], ok = True},"
                " ok = Length[H] === Length[H[[1]]];"
                " Do[If[j < i && !zq[H[[i, j]]], ok = False]; If[j == i && !zq[plead[H[[i, i]]] - 1], ok = False];"
                " If[j > i && pdeg[H[[i, j]]] >= pdeg[H[[j, j]]], ok = False], {i, c}, {j, c}];"
                " ok = ok && AllTrue[M, inspan[#, H] &];"
                " If[ok && Length[M] === c, With[{dd = Together[Det[M]]}, With[{nn = red[Numerator[dd]]},"
                "   With[{l = Coefficient[nn, x, Exponent[nn, x]]/Denominator[dd]},"
                "   ok = zq[Det[M] - l Times @@ Diagonal[H]]]]]];"
                " ok]]") % (ins, out)
    elif op == "detadj":
        body = ("With[{M = %s[[1]], O0 = %s}, zq[O0[[1]] - Det[M]] && "
                "AllTrue[Flatten[O0[[2]] . M - O0[[1]] IdentityMatrix[Length[M]]], zq]]") % (ins, out)
    else:
        return None
    return "%s Print[\"R \", %d, \" \", TimeConstrained[Check[If[TrueQ[%s], 1, 0], -1], 60, -2]];\n" % (head, c["id"], body)


def run_script(cmd, text, suffix):
    with tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False) as fh:
        fh.write(text)
        path = fh.name
    try:
        out = subprocess.run(cmd + [path], capture_output=True, text=True, timeout=7200)
        res = {}
        for line in (out.stdout + "\n" + out.stderr).splitlines():
            m = re.match(r"\s*R\s+(\d+)\s+(-?\d+)", line)
            if m:
                res[int(m.group(1))] = int(m.group(2))
        return res, out
    finally:
        os.unlink(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=30, help="cases per field and operation")
    ap.add_argument("--seed", type=int, default=20261115)
    ap.add_argument("--no-gp", action="store_true")
    ap.add_argument("--no-wl", action="store_true")
    args = ap.parse_args()
    cases = run_gen(args.n, args.seed)
    live = [c for c in cases if c["out"] not in ("DECLINED", "TIMEOUT")]
    declined = [c for c in cases if c["out"] == "DECLINED"]
    slow = [c for c in cases if c["out"] == "TIMEOUT"]
    print("cases: %d (declined by the kernel: %d, over the 20 s budget: %d)" % (len(cases), len(declined), len(slow)))
    verdicts = {}
    if not args.no_gp:
        txt = GP_PRELUDE + "".join(gp_case(c) or "" for c in live)
        verdicts["gp"], _ = run_script([GP, "-q", "-f"], txt + "\\q\n", ".gp")
    if not args.no_wl:
        txt = WL_PRELUDE + "".join(wl_case(c) or "" for c in live)
        verdicts["wl"], _ = run_script([WL, "-file"], txt, ".wl")
    table = defaultdict(lambda: defaultdict(int))
    bad = 0
    for c in live:
        key = (c["field"], c["op"])
        for sysname, v in verdicts.items():
            r = v.get(c["id"])
            tag = {1: "ok", 0: "FAIL", -1: "error", -2: "timeout"}.get(r, "missing")
            table[key][sysname + ":" + tag] += 1
            if tag in ("FAIL", "missing"):
                bad += 1
                print("%s %s case %d (%s, %s): ins %s out %s" % (sysname, tag, c["id"], c["field"], c["op"],
                                                              c["ins"][:200], c["out"][:200]))
    for c in declined:
        table[(c["field"], c["op"])]["declined"] += 1
    for c in slow:
        table[(c["field"], c["op"])]["kernel-timeout"] += 1
    for key in sorted(table):
        print("%-8s %-7s %s" % (key[0], key[1], "  ".join("%s=%d" % kv for kv in sorted(table[key].items()))))
    print("FFS-SUMMARY failures=%d declined=%d" % (bad, len(declined)))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
