#!/usr/bin/env python3
"""tools/ffield_bench/bench.py -- number-field polynomial arithmetic, timed in
Mathilda, PARI/GP and the reference CAS (a developer tool, not a CI gate: the
reference CAS is not available in CI).

Fields K = Q[t]/(mu), each with a radical spelling of t for the user-level
calls; for each field and size, random polynomials over K in x with a planted
common factor.  Operations:

  mul     A*B                           (user level: Expand over the radicals)
  divrem  A = q B + r
  gcd     the monic gcd                 (planted factor of degree deg/3)
  xgcd    g, s, t with s A + t B = g
  norm    Resultant[A, mu, t]           (the norm of A down to Q[x])
  factor  factorisation over K          (A = planted factor times cofactor)

Two levels in Mathilda: the kernel (FLINT` builtins with Extension -> {t, mu},
coefficients as polynomials in t) and the user level (PolynomialGCD / Factor /
Expand on the radical spelling).  PARI works on POLMOD coefficients; the
reference CAS on the radical spelling (it has no xgcd over algebraic
coefficients).  Each system times only the operation, on inputs parsed
beforehand, repeating it until 0.2 s have passed.  A timeout or an unevaluated
result is reported as such, and every system's gcd degree is reported so the
answers can be cross-checked.

  python3 tools/ffield_bench/bench.py [--quick] [--out results.json]
"""
import argparse, glob, json, os, random, shutil, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MATH = os.path.join(ROOT, "Mathilda")

FIELDS = [
    # name, mu (in t), radical spelling of t
    ("sqrt2", "t^2 - 2", "Sqrt[2]"),
    ("omega", "t^2 + t + 1", "(-1 + I*Sqrt[3])/2"),
    ("cbrt2", "t^3 - 2", "2^(1/3)"),
    ("sqrt2sqrt3", "t^4 - 10*t^2 + 1", "(Sqrt[2] + Sqrt[3])"),
]
SIZES = [(8, 2), (24, 2), (24, 20), (48, 4)]   # (degree in x, coefficient digits)
QUICK_SIZES = [(8, 2), (24, 20)]
OPS = ["mul", "divrem", "gcd", "xgcd", "norm", "factor"]


def find_wl():
    env = os.environ.get("FFS_WL")
    if env and os.path.exists(env):
        return env
    w = shutil.which("wolframscript")
    if w:
        return w
    c = sorted(glob.glob("/Applications/*.app/Contents/MacOS/wolframscript"))
    return c[0] if c else None


def rand_poly(rng, deg, n, digits):
    """a random polynomial over K: sum_d (sum_j c_dj t^j) x^d, monic-ish in x"""
    terms = []
    for d in range(deg + 1):
        for j in range(n):
            c = rng.randint(-10 ** digits, 10 ** digits)
            if c:
                terms.append("(%d)*t^%d*x^%d" % (c, j, d))
    terms.append("x^%d" % (deg + 1))
    return " + ".join(terms)


def make_cases(sizes, seed):
    rng = random.Random(seed)
    cases = []
    for fname, mu, rad in FIELDS:
        n = int(mu.split("^")[1].split()[0])
        for deg, dig in sizes:
            k = max(1, deg // 3)
            C = rand_poly(rng, k, n, dig)
            P = rand_poly(rng, deg - k, n, dig)
            Q = rand_poly(rng, deg - k, n, dig)
            cases.append({"id": "%s-d%d-c%d" % (fname, deg, dig), "field": fname, "mu": mu,
                          "rad": rad, "n": n, "deg": deg, "digits": dig,
                          "C": C, "P": P, "Q": Q})
    return cases


def mathilda_script(cases, outp):
    L = ['tm[f_] := Module[{t0 = AbsoluteTime[], k = 0, r}, While[k == 0 || AbsoluteTime[] - t0 < 0.2, r = f[]; k++]; {(AbsoluteTime[] - t0)/k, r}];',
         'tc[f_] := TimeConstrained[tm[f], 60, {"TIMEOUT", $Failed}];',
         'red[e_, mu_] := Expand[PolynomialRemainder[Expand[e], mu, t]];',
         'row[id_, op_, v_, info_] := Print["ROW\t", id, "\t", op, "\t", If[NumberQ[v[[1]]], ToString[CForm[N[v[[1]]]]], v[[1]]], "\t", If[v[[1]] === "TIMEOUT", "-", Quiet[ToString[info[v[[2]]]]]]];']
    for c in cases:
        cid = c["id"]
        L.append('Module[{mu = %s, rad = %s, C, P, Q, A, B, Ar, Br},' % (c["mu"], c["rad"]))
        L.append('  C = %s; P = %s; Q = %s;' % (c["C"], c["P"], c["Q"]))
        L.append('  A = red[C P, mu]; B = red[C Q, mu]; Ar = Expand[A /. t -> rad]; Br = Expand[B /. t -> rad];')
        L.append('  row["%s", "mul", tc[Function[Expand[Ar Br]]], 0 &];' % cid)
        L.append('  row["%s", "divrem", tc[Function[FLINT`PolynomialQuotientRemainder[A, B, x, Extension -> {t, mu}]]], 0 &];' % cid)
        L.append('  row["%s", "gcd", tc[Function[FLINT`PolynomialFieldGCD[A, B, x, Extension -> {t, mu}]]], Exponent[#, x] &];' % cid)
        L.append('  row["%s", "gcd_user", tc[Function[PolynomialGCD[Ar, Br, Extension -> Automatic]]], Exponent[#, x] &];' % cid)
        L.append('  row["%s", "xgcd", tc[Function[FLINT`PolynomialExtendedGCD[A, B, x, Extension -> {t, mu}]]], Exponent[#[[1]], x] &];' % cid)
        L.append('  row["%s", "norm", tc[Function[Resultant[A, mu, t]]], 0 &];' % cid)
        L.append('  row["%s", "factor_user", tc[Function[Factor[Ar, Extension -> rad]]], Length[FactorList[#, Extension -> rad]] &]];' % cid)
    return "\n".join(L)


def gp_script(cases, outp):
    L = ['default(parisize, 1000000000); default(realprecision, 10);',
         'tm(f) = my(t0 = getabstime(), k = 0, r); while(k == 0 || getabstime() - t0 < 200, r = f(); k++); [(getabstime() - t0)/1000./k, r];',
         'out = [];']
    for c in cases:
        mu = c["mu"].replace("t", "y")
        conv = lambda s: s.replace("t", "y")
        L.append('{ my(mu = %s, C = %s, P = %s, Q = %s, A, B, nf, v1, v2, v3, v4, v5, v6);' %
                 (mu, conv(c["C"]), conv(c["P"]), conv(c["Q"])))
        L.append('  A = C*P*Mod(1, mu); B = C*Q*Mod(1, mu); nf = nfinit(mu);')
        L.append('  v1 = tm(()->A*B); v2 = tm(()->divrem(A, B)); v3 = tm(()->gcd(A, B)); v4 = tm(()->gcdext(A, B));')
        L.append('  v5 = tm(()->polresultant(lift(A), mu, y)); v6 = tm(()->nffactor(nf, lift(A)));')
        cid = c["id"]
        L.append('  print("ROW\\t%s\\tmul\\t", v1[1], "\\t0"); print("ROW\\t%s\\tdivrem\\t", v2[1], "\\t0");' % (cid, cid))
        L.append('  print("ROW\\t%s\\tgcd\\t", v3[1], "\\t", poldegree(v3[2])); print("ROW\\t%s\\txgcd\\t", v4[1], "\\t", poldegree(v4[2][3]));' % (cid, cid))
        L.append('  print("ROW\\t%s\\tnorm\\t", v5[1], "\\t0"); print("ROW\\t%s\\tfactor\\t", v6[1], "\\t", #v6[2][,1]); }' % (cid, cid))
    L.append('quit;')
    return "\n".join(L)


def wl_script(cases, outp):
    L = ['tm[f_] := Module[{t0 = AbsoluteTime[], k = 0, r}, While[k == 0 || AbsoluteTime[] - t0 < 0.2, r = f[]; k++]; {(AbsoluteTime[] - t0)/k, r}];',
         'tc[f_] := TimeConstrained[tm[f], 60, {"TIMEOUT", $Failed}];',
         'red[e_, mu_] := Expand[PolynomialRemainder[Expand[e], mu, t]];',
         'row[id_, op_, v_, info_] := Print["ROW\t", id, "\t", op, "\t", If[NumberQ[v[[1]]], ToString[CForm[N[v[[1]]]]], v[[1]]], "\t", If[v[[1]] === "TIMEOUT", "-", Quiet[ToString[info[v[[2]]]]]]];']
    for c in cases:
        L.append('Module[{mu = %s, rad = %s, C, P, Q, A, B, Ar, Br},' % (c["mu"], c["rad"]))
        L.append('  C = %s; P = %s; Q = %s;' % (c["C"], c["P"], c["Q"]))
        L.append('  A = red[C P, mu]; B = red[C Q, mu]; Ar = Expand[A /. t -> rad]; Br = Expand[B /. t -> rad];')
        cid = c["id"]
        L.append('  row["%s", "mul", tc[Function[Expand[Ar Br]]], 0 &];' % cid)
        L.append('  row["%s", "divrem", tc[Function[PolynomialQuotientRemainder[Ar, Br, x]]], 0 &];' % cid)
        L.append('  row["%s", "gcd_user", tc[Function[PolynomialGCD[Ar, Br, Extension -> Automatic]]], Exponent[#, x] &];' % cid)
        L.append('  row["%s", "norm", tc[Function[Resultant[A, mu, t]]], 0 &];' % cid)
        L.append('  row["%s", "factor_user", tc[Function[Factor[Ar, Extension -> rad]]], Length[FactorList[#, Extension -> rad]] &]];' % cid)
    return "\n".join(L)


def run(cmd, timeout):
    try:
        p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=timeout,
                           stdin=subprocess.DEVNULL)
        return p.returncode, p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        return -1, "TIMEOUT"


def fmt(v):
    if v is None:
        return "-"
    if isinstance(v, str):
        return v
    if isinstance(v, (int, float)):
        if v >= 1:
            return "%.2f s" % v
        return "%.2f ms" % (1000 * v)
    return str(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=20261115)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    cases = make_cases(QUICK_SIZES if args.quick else SIZES, args.seed)
    tmp = tempfile.mkdtemp(prefix="ffbench")
    files = {k: os.path.join(tmp, k) for k in ["m.m", "g.gp", "w.wls", "m.json", "g.txt", "w.json"]}
    open(files["m.m"], "w").write(mathilda_script(cases, files["m.json"]))
    open(files["g.gp"], "w").write(gp_script(cases, files["g.txt"]))
    open(files["w.wls"], "w").write(wl_script(cases, files["w.json"]))
    results = {}
    print("scripts in", tmp, file=sys.stderr)

    def rows(out):
        r = {}
        for line in out.splitlines():
            if not line.startswith("ROW\t"):
                continue
            f = line.split("\t")
            if len(f) < 5:
                continue
            cid, op, tm, info = f[1].strip(), f[2].strip(), f[3].replace(" ", ""), f[4].strip()
            try:
                tv = float(tm)
            except ValueError:
                tv = tm
            r.setdefault(cid, {})[op] = [tv, info]
        return r

    rc, out = run([MATH, "-file", files["m.m"]], 7200)
    results["mathilda"] = rows(out)
    if not results["mathilda"]:
        print("Mathilda run failed:\n" + out[-2000:], file=sys.stderr)
    if shutil.which("gp"):
        rc, out = run(["gp", "-q", files["g.gp"]], 7200)
        results["pari"] = rows(out)
        if not results["pari"]:
            print("PARI run failed:\n" + out[-2000:], file=sys.stderr)
    wl = find_wl()
    if wl:
        rc, out = run([wl, "-file", files["w.wls"]], 7200)
        results["refcas"] = rows(out)
    # table
    print("%-22s %-8s %12s %12s %12s %12s" % ("case", "op", "Mathilda", "M (user)", "PARI", "ref CAS"))
    for c in cases:
        cid = c["id"]
        m = results.get("mathilda", {}).get(cid, {})
        g = results.get("pari", {}).get(cid, {})
        w = results.get("refcas", {}).get(cid, {})
        def first(x):
            return x[0] if isinstance(x, list) else x
        table = [("mul", None, first(m.get("mul")), first(g.get("mul")), first(w.get("mul"))),
                 ("divrem", first(m.get("divrem")), None, first(g.get("divrem")), first(w.get("divrem"))),
                 ("gcd", first(m.get("gcd")), first(m.get("gcd_user")), first(g.get("gcd")), first(w.get("gcd_user"))),
                 ("xgcd", first(m.get("xgcd")), None, first(g.get("xgcd")), None),
                 ("norm", first(m.get("norm")), None, first(g.get("norm")), first(w.get("norm"))),
                 ("factor", None, first(m.get("factor_user")), first(g.get("factor")), first(w.get("factor_user")))]
        for op, a, b, pg, ww in table:
            print("%-22s %-8s %12s %12s %12s %12s" % (cid, op, fmt(a), fmt(b), fmt(pg), fmt(ww)))
        info = lambda d, k: (d.get(k) or [None, None])[1]
        print("%-22s gcd deg (M, M user, PARI, ref): %s %s %s %s | xgcd deg M/PARI %s %s | factors M/PARI/ref %s %s %s" % (
            "", info(m, "gcd"), info(m, "gcd_user"), info(g, "gcd"), info(w, "gcd_user"),
            info(m, "xgcd"), info(g, "xgcd"), info(m, "factor_user"), info(g, "factor"), info(w, "factor_user")))
    if args.out:
        json.dump({"cases": [{k: v for k, v in c.items() if k not in ("C", "P", "Q")} for c in cases],
                   "results": results}, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
