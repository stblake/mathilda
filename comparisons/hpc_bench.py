#!/usr/bin/env python3
"""
Mathilda vs Mathematica on classical high-performance-computing kernels.

Runs each benchmark in both systems, in one process per system, timing with
each system's AbsoluteTiming (elapsed wall clock -- NOT Timing[], which reports
CPU time summed over threads and so over-reports every threaded NDArray path and
every BLAS call by roughly the core count).

Both systems run the SAME source text for every benchmark; a `wl` override is
only used where the two languages genuinely differ. Where a benchmark has a
cheap scalar answer, that answer is recorded from both systems and compared --
a timing row is only meaningful if the two agree.

Usage:
    comparisons/hpc_bench.py                 # full run
    comparisons/hpc_bench.py --scale 0.1     # smaller sizes, for a smoke test
    comparisons/hpc_bench.py --only fft,sort # a subset, by id
    comparisons/hpc_bench.py --system mathilda
Output: a markdown table on stdout, and the raw JSON to --json.
"""
import argparse, json, os, re, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MATHILDA = os.environ.get("MATHILDA_BIN", os.path.join(ROOT, "Mathilda"))
WOLFRAM = os.environ.get(
    "WOLFRAMSCRIPT", "/Applications/Mathematica.app/Contents/MacOS/wolframscript")

# --------------------------------------------------------------------------
# Benchmarks.
#
# id      short key for --only
# group   section heading in the report
# name    row label
# setup   run once, untimed (data construction, kernel definition)
# expr    the timed expression; ends in ';' where the value is a big array
# check   optional: a cheap scalar recorded from both systems and compared
# reps    timed repetitions; the reported number is the MINIMUM
# scale   which size knobs to shrink under --scale
# --------------------------------------------------------------------------
B = []


def bench(**kw):
    kw.setdefault("reps", 5)
    kw.setdefault("check", None)
    kw.setdefault("wl", None)
    kw.setdefault("wl_setup", None)
    kw.setdefault("note", "")
    kw.setdefault("cold", False)
    B.append(kw)


LA = "Dense linear algebra (BLAS / LAPACK)"
SP = "Spectral"
AR = "Array and memory primitives (n = 10^7)"
ST = "Stencils and PDE relaxation"
SC = "Scalar kernels via Compile[]"
IN = "Integer and combinatorial"
BN = "Arbitrary precision (GMP / MPFR)"

# ---- dense linear algebra -------------------------------------------------
bench(id="matmul", group=LA, name="Matrix multiply, 1000x1000",
      setup="n=1000; A=RandomReal[{0,1},{n,n}]; Bm=RandomReal[{0,1},{n,n}];",
      expr="A . Bm;", reps=5)
bench(id="solve", group=LA, name="LinearSolve, 1000x1000",
      setup="n=1000; A=RandomReal[{0,1},{n,n}]; bv=RandomReal[{0,1},n];",
      expr="LinearSolve[A,bv];", reps=5)
bench(id="inverse", group=LA, name="Inverse, 500x500",
      setup="n=500; A=RandomReal[{0,1},{n,n}];", expr="Inverse[A];", reps=5)
bench(id="det", group=LA, name="Det, 500x500",
      setup="n=500; A=RandomReal[{0,1},{n,n}];", expr="Det[A];", reps=5)
bench(id="qr", group=LA, name="QRDecomposition, 500x500",
      setup="n=500; A=RandomReal[{0,1},{n,n}];", expr="QRDecomposition[A];", reps=3)
bench(id="eigen", group=LA, name="Eigenvalues, 300x300 symmetric",
      setup="n=300; A=RandomReal[{0,1},{n,n}]; S=(A+Transpose[A])/2.;",
      expr="Eigenvalues[S];", reps=3)
bench(id="svd", group=LA, name="SingularValueDecomposition, 300x300",
      setup="n=300; A=RandomReal[{0,1},{n,n}];",
      expr="SingularValueDecomposition[A];", reps=3)

# ---- spectral -------------------------------------------------------------
bench(id="fft", group=SP, name="Fourier, 2^20 reals",
      setup="v=RandomReal[{0,1},2^20];", expr="Fourier[v];", reps=5)
bench(id="fft2", group=SP, name="Fourier, 2^20 complex",
      setup="vz=RandomReal[{0,1},2^20] + I RandomReal[{0,1},2^20];",
      expr="Fourier[vz];", reps=5)

# ---- array / memory primitives -------------------------------------------
ARSET = "n=10^7; x=RandomReal[{0,1},n]; y=RandomReal[{0,1},n];"
bench(id="triad", group=AR, name="STREAM triad, a = b + 3 c",
      setup=ARSET, expr="x + 3. y;", reps=5)
bench(id="total", group=AR, name="Total (reduction)", setup=ARSET,
      expr="Total[x];", reps=5)
bench(id="accum", group=AR, name="Accumulate (prefix scan)", setup=ARSET,
      expr="Accumulate[x];", reps=5)
bench(id="sort", group=AR, name="Sort", setup=ARSET, expr="Sort[x];", reps=3)
bench(id="sin", group=AR, name="Sin (elementwise)", setup=ARSET,
      expr="Sin[x];", reps=5)
bench(id="exp", group=AR, name="Exp (elementwise)", setup=ARSET,
      expr="Exp[x];", reps=5)
bench(id="dot", group=AR, name="Dot (inner product)", setup=ARSET,
      expr="x . y;", reps=5)
# The structural family. Every one of these was on the marshalling path (one Expr
# per element) until this round of work; they are the operations array-style
# numerical code is actually written out of.
bench(id="rotate", group=AR, name="RotateLeft", setup=ARSET,
      expr="RotateLeft[x,3];", reps=5)
bench(id="reverse", group=AR, name="Reverse", setup=ARSET,
      expr="Reverse[x];", reps=5)
bench(id="joinv", group=AR, name="Join (two 10^7 vectors)", setup=ARSET,
      expr="Join[x,y];", reps=3)
bench(id="part2", group=AR, name="Partition, length 2", setup=ARSET,
      expr="Partition[x,2];", reps=3)
bench(id="diffs", group=AR, name="Differences", setup=ARSET,
      expr="Differences[x];", reps=3)
bench(id="riffle", group=AR, name="Riffle with a Real", setup=ARSET,
      expr="Riffle[x,0.];", reps=3)
bench(id="padr", group=AR, name="PadRight, Real fill", setup=ARSET,
      expr="PadRight[x,n+1,0.];", reps=3)
bench(id="padrd", group=AR, name="PadRight, default (exact 0) fill", setup=ARSET,
      expr="PadRight[x,n+1];", reps=3,
      note="mixed exact/inexact result: unpacked in BOTH systems, by definition")

# ---- stencils -------------------------------------------------------------
bench(id="jacobi", group=ST, name="Jacobi 5-point relaxation, 512^2, 100 sweeps",
      setup="n=512; u0=RandomReal[{0,1},{n,n}]; "
            "jac[u_]:=(RotateLeft[u,{1,0}]+RotateRight[u,{1,0}]"
            "+RotateLeft[u,{0,1}]+RotateRight[u,{0,1}])/4.;",
      expr="Nest[jac,u0,100];", reps=3)
bench(id="life", group=ST, name="Game of Life, 256^2, 100 generations",
      setup="n=256; g0=RandomInteger[1,{n,n}]; "
            "nb[g_]:=Sum[RotateLeft[g,{i,j}],{i,-1,1},{j,-1,1}]-g; "
            "life[g_]:=With[{k=nb[g]}, "
            "UnitStep[k-3]UnitStep[3-k] + UnitStep[k-2]UnitStep[2-k] g];",
      expr="Nest[life,g0,100];", reps=3)

# ---- scalar kernels via Compile[] ----------------------------------------
bench(id="logistic", group=SC, name="Logistic map, 10^7 iterations",
      setup="lg=Compile[{{x0,_Real},{m,_Integer}}, "
            "Module[{x=x0,k=0}, While[k<m, x=3.9 x (1.-x); k=k+1]; x]];",
      expr="lg[0.5, 10^7];", check="lg[0.5, 10^7]", reps=3)
bench(id="mandel", group=SC, name="Mandelbrot, 800x800, 100 iterations",
      setup="mandel=Compile[{{cx,_Real},{cy,_Real},{mx,_Integer}}, "
            "Module[{zx=0.,zy=0.,t=0.,k=0}, "
            "While[k<mx && zx zx + zy zy < 4., "
            "t=zx zx - zy zy + cx; zy=2. zx zy + cy; zx=t; k=k+1]; k]]; "
            "n=800; h=2.5/(n-1);",
      expr="Table[mandel[x,yy,100],{yy,-1.25,1.25,h},{x,-2.,0.5,h}];",
      check="Total[Table[mandel[x,yy,100],{yy,-1.25,1.25,h},{x,-2.,0.5,h}],2]",
      reps=3)
bench(id="mcpi", group=SC, name="Monte Carlo pi, 10^7 samples (vectorized)",
      setup="n=10^7;",
      expr="4. Total[UnitStep[1. - RandomReal[{0,1},n]^2 "
           "- RandomReal[{0,1},n]^2]]/n;", reps=3)
# The point set is a cubic LATTICE, deliberately. A pseudo-random one built from
# Mod[a k, m] cycles, so with 3n coordinates it repeats POINTS; two coincident
# bodies make r = 0, the 1/r^12 term overflows to Inf, and the compiled call then
# correctly bails to the interpreter on the non-finite result -- which read as an
# 89x performance cliff between n=1000 and n=1200 and answered Indeterminate.
# A lattice with spacing 1.2 has no coincident pairs and a finite energy.
bench(id="lj", group=SC, name="Lennard-Jones energy, 1452 bodies (all pairs)",
      setup="nb=1452; pts=Flatten[Table[1.2 {i,j,k2}, "
            "{i,0,10},{j,0,10},{k2,0,11}]]; "
            "lj=Compile[{{p,_Real,1},{m,_Integer}}, "
            "Module[{e=0.,i=0,j=0,dx=0.,dy=0.,dz=0.,r2=0.,ir6=0.}, "
            "For[i=0,i<m-1,i++, For[j=i+1,j<m,j++, "
            "dx=p[[3i+1]]-p[[3j+1]]; dy=p[[3i+2]]-p[[3j+2]]; "
            "dz=p[[3i+3]]-p[[3j+3]]; r2=dx dx+dy dy+dz dz; "
            "ir6=1./(r2 r2 r2); e=e+4. ir6 (ir6-1.)]]; e]];",
      expr="lj[pts, nb];", check="lj[pts, nb]", reps=3)

# ---- integer / combinatorial --------------------------------------------
bench(id="sieve", group=IN, name="Sieve of Eratosthenes to 10^7",
      setup="sv=Compile[{{m,_Integer}}, Module[{s=ConstantArray[1,m],i=2,j=0,c=0}, "
            "While[i i<=m, If[s[[i]]==1, j=i i; While[j<=m, s[[j]]=0; j=j+i]]; "
            "i=i+1]; i=2; While[i<=m, c=c+s[[i]]; i=i+1]; c]];",
      expr="sv[10^7];", check="sv[10^7]", reps=3)
bench(id="collatz", group=IN, name="Collatz longest chain below 10^6",
      setup="cz=Compile[{{m,_Integer}}, Module[{bl=0,k=1,q=0,len=0}, "
            "While[k<=m, q=k; len=0; While[q>1, "
            "If[Mod[q,2]==0, q=Quotient[q,2], q=3q+1]; len=len+1]; "
            "If[len>bl, bl=len]; k=k+1]; bl]];",
      expr="cz[10^6];", check="cz[10^6]", reps=3)
bench(id="fib", group=IN, name="Naive recursive Fibonacci, fib(25)",
      setup="Clear[fib]; fib[0]=0; fib[1]=1; fib[k_]:=fib[k-1]+fib[k-2];",
      expr="fib[25];", check="fib[25]", reps=3,
      note="rule dispatch, not compiled in either system")
bench(id="primepi", group=IN, name="PrimePi[10^9]", setup="",
      expr="PrimePi[10^9 + kk];", check="PrimePi[10^9]", reps=3)

# ---- arbitrary precision ------------------------------------------------
bench(id="pi", group=BN, name="pi to 100,000 digits", setup="",
      expr="N[Pi,100000 + kk];", check="Round[10^20 (N[Pi,100000] - 3)]", reps=3,
      cold=True)
bench(id="fact", group=BN, name="50000! (exact)", setup="",
      expr="(50000 + kk)!;", check="IntegerLength[50000!]", reps=3)
# IntegerLength[], not `expr;`. Wolfram does not materialise a discarded bignum
# product: `(p + kk) q;` reads 138 us and the forced form 2.96 ms, which would
# have been reported as a 24x gap that does not exist. Array ops were checked in
# both forms and agree within noise; bignum arithmetic is the exception.
bench(id="bigmul", group=BN, name="Product of two 10^6-bit integers",
      setup="p=2^1000003 - 1; q=3^631305;",
      expr="IntegerLength[(p + kk) q];", check="IntegerLength[p q]", reps=5)
bench(id="bigpow", group=BN, name="3^1000000", setup="",
      expr="3^(1000000 + kk);", check="IntegerLength[3^1000000]", reps=3)

SCALED = {  # id -> (pattern, replacement-fn) applied to setup+expr under --scale
    "matmul": [("n=1000", "n")], "solve": [("n=1000", "n")],
    "inverse": [("n=500", "n")], "det": [("n=500", "n")],
    "qr": [("n=500", "n")], "eigen": [("n=300", "n")], "svd": [("n=300", "n")],
}


def apply_scale(b, s):
    """Shrink the dominant size knob so the whole suite can smoke-test fast."""
    if s == 1.0:
        return b
    b = dict(b)
    for key in ("setup", "expr", "check", "wl", "wl_setup"):
        t = b.get(key)
        if not t:
            continue
        def shrink(m):
            v = int(m.group(2))
            return m.group(1) + str(max(2, int(v * (s ** 0.5) if v > 64 else v)))
        t = re.sub(r"(n=|nb=|m=)(\d+)", shrink, t)
        t = t.replace("10^7", "10^5").replace("10^6", "10^4")
        t = t.replace("2^20", "2^14").replace("100000", "2000")
        t = t.replace("50000!", "2000!").replace("3^1000000", "3^20000")
        t = t.replace("10^9", "10^6").replace(",100]", ",20]")
        t = t.replace("u0,100]", "u0,5]").replace("g0,100]", "g0,5]")
        b[key] = t
    return b


# --------------------------------------------------------------------------
# Mathilda: NDJSON pipe protocol (non-TTY stdin).
# --------------------------------------------------------------------------
def run_mathilda(benches):
    reqs, tags = [], []

    def send(expr, tag=None):
        reqs.append({"id": len(reqs), "expr": expr})
        tags.append(tag)

    for b in benches:
        if b["setup"]:
            send(b["setup"])
        # warm-up, untimed: first touch pays page faults and any one-off compile.
        # Skipped for a `cold` benchmark: some builtins memoise so hard that no
        # amount of argument variation defeats them (Mathilda's Pi cache keeps
        # guard digits, so asking for 100003 after 100000 is a cache hit and the
        # timed reps read 1 us). For those the FIRST, cold call is the honest
        # measurement of the work, so it is the one reported.
        if not b["cold"]:
            send("kk = 0; " + b["expr"])
        # {min, max} over reps. The MAX is the caching tripwire: a builtin with an
        # internal cache (PrimePi's table, Pi's digits) answers the 2nd rep in
        # ~0 s, and reporting the min alone would make the row vacuous. `kk` is
        # bound so an expression can vary its argument per rep and defeat that.
        send("tsq = Table[First[AbsoluteTiming[%s]],{kk,1,%d}]; {Min[tsq], Max[tsq]}"
             % (b["expr"], 1 if b["cold"] else b["reps"]), ("t", b["id"]))
        if b["check"]:
            send("InputForm[%s]" % b["check"], ("v", b["id"]))
    send("Quit[]")

    payload = "".join(json.dumps(r) + "\n" for r in reqs) + \
              json.dumps({"type": "quit"}) + "\n"
    t0 = time.time()
    p = subprocess.run([MATHILDA], input=payload, capture_output=True, text=True)
    out = {}
    for ln in p.stdout.splitlines():
        ln = ln.strip()
        if not ln.startswith("{"):
            continue
        try:
            o = json.loads(ln)
        except Exception:
            continue
        if o.get("type") != "expr":
            continue
        tag = tags[o["id"]] if o["id"] < len(tags) else None
        if tag:
            out.setdefault(tag[1], {})[tag[0]] = o.get("payload", "").strip()
    sys.stderr.write("[mathilda] %.1fs wall\n" % (time.time() - t0))
    if p.returncode != 0:
        sys.stderr.write("[mathilda] exit %d\n%s\n" % (p.returncode, p.stderr[-2000:]))
    return out


# --------------------------------------------------------------------------
# Mathematica: one wolframscript file (startup is ~3 s, so batch everything).
# --------------------------------------------------------------------------
def run_wolfram(benches, tmpdir):
    lines = ["$HistoryLength = 0;"]
    for b in benches:
        setup = b["wl_setup"] or b["setup"]
        expr = b["wl"] or b["expr"]
        if setup:
            lines.append(setup)
        if not b["cold"]:
            lines.append("kk = 0; " + expr)                    # warm-up
        lines.append('tsq = Table[First[AbsoluteTiming[%s]],{kk,1,%d}]; '
                     'Print["@t|%s|", {Min[tsq], Max[tsq]}];'
                     % (expr, 1 if b["cold"] else b["reps"], b["id"]))
        if b["check"]:
            chk = b["check"] if not b["wl"] else b["check"]
            lines.append('Print["@v|%s|", InputForm[%s]];' % (b["id"], chk))
    path = os.path.join(tmpdir, "hpc_bench.wls")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    t0 = time.time()
    p = subprocess.run([WOLFRAM, "-file", path], capture_output=True, text=True)
    out = {}
    for ln in p.stdout.splitlines():
        m = re.match(r"@([tv])\|([^|]+)\|(.*)", ln.strip())
        if m:
            v = m.group(3).strip()
            # Print[InputForm[x]] emits the literal wrapper "InputForm[...]" in
            # wolframscript, so every value comparison read as a mismatch. Strip it.
            if v.startswith("InputForm[") and v.endswith("]"):
                v = v[len("InputForm["):-1]
            out.setdefault(m.group(2), {})[m.group(1)] = v
    sys.stderr.write("[wolfram] %.1fs wall\n" % (time.time() - t0))
    if p.stderr.strip():
        sys.stderr.write("[wolfram] stderr:\n" + p.stderr[-2000:] + "\n")
    return out


def values_agree(a, b):
    """Do two printed values agree?

    Exact match first, then a numeric comparison to the SHORTER of the two
    printings: Mathilda's printer emits 6 significant figures where Wolfram emits
    17, so -5317.09 and -5317.087608456237 are the same answer and a string
    comparison would call them different.
    """
    if a == b:
        return True
    try:
        x, y = float(a.replace("*^", "e")), float(b.replace("*^", "e"))
    except ValueError:
        return False
    if x == y:
        return True
    # Significant digits of the SHORTER printing: strip sign, decimal point and
    # leading zeros (0.123564 has six significant digits, not seven), then allow
    # half a unit in the last of them, since printing rounds.
    def sigdigits(t):
        d = t.replace("-", "").replace(".", "").lstrip("0")
        return len(d.rstrip("0")) or 1
    sig = min(sigdigits(a), sigdigits(b))
    tol = 5.0 * 10.0 ** -sig
    return abs(x - y) <= tol * max(abs(x), abs(y), 1e-300)


def as_pair(s):
    """Parse a printed {min, max} into a (min, max) float pair."""
    if not s:
        return (None, None)
    t = s.replace("*^", "e").replace("`", "").strip().strip("{}")
    parts = [x.strip() for x in t.split(",")]
    try:
        v = [float(x) for x in parts]
    except ValueError:
        return (None, None)
    if len(v) == 1:
        return (v[0], v[0])
    return (v[0], v[1])


def fmt(t):
    if t is None:
        return "—"
    if t < 1e-3:
        return "%.0f us" % (t * 1e6)
    if t < 1.0:
        return "%.2f ms" % (t * 1e3)
    return "%.3f s" % t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--only", default="")
    ap.add_argument("--system", default="both",
                    choices=["both", "mathilda", "wolfram"])
    ap.add_argument("--json", default="")
    a = ap.parse_args()

    sel = set(x.strip() for x in a.only.split(",") if x.strip())
    benches = [apply_scale(b, a.scale) for b in B if not sel or b["id"] in sel]
    tmp = os.environ.get("CLAUDE_JOB_DIR", "/tmp")
    tmp = os.path.join(tmp, "tmp") if os.path.isdir(os.path.join(tmp, "tmp")) else tmp

    md = run_mathilda(benches) if a.system in ("both", "mathilda") else {}
    wl = run_wolfram(benches, tmp) if a.system in ("both", "wolfram") else {}

    rows, group = [], None
    print("| Benchmark | Mathilda | Mathematica 14.0 | ratio |")
    print("|---|---:|---:|---:|")
    for b in benches:
        i = b["id"]
        (tm, tmx) = as_pair(md.get(i, {}).get("t"))
        (tw, twx) = as_pair(wl.get(i, {}).get("t"))
        vm, vw = md.get(i, {}).get("v"), wl.get(i, {}).get("v")
        if b["group"] != group:
            group = b["group"]
            print("| **%s** | | | |" % group)
        ratio = "—"
        if tm and tw:
            ratio = ("%.2fx" % (tw / tm)) if tw >= tm else ("1/%.2fx" % (tm / tw))
        flag = ""
        if b["check"] and vm and vw and not values_agree(vm, vw):
            flag = "  ⚠ VALUE MISMATCH"
        # min far below max => an internal cache answered the later reps
        for who, lo, hi in (("mathilda", tm, tmx), ("wolfram", tw, twx)):
            if lo and hi and lo > 0 and hi / lo > 5:
                flag += f"  ⚠ {who} max/min={hi/lo:.0f}x (cache?)"
        print("| %s | %s | %s | %s |%s" % (b["name"], fmt(tm), fmt(tw), ratio, flag))
        rows.append(dict(id=i, name=b["name"], group=b["group"],
                         mathilda=tm, mathilda_max=tmx,
                         wolfram=tw, wolfram_max=twx, vm=vm, vw=vw,
                         agree=(None if not b["check"] else (vm == vw))))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(rows, f, indent=1)
    bad = [r for r in rows if r["agree"] is False]
    if bad:
        sys.stderr.write("\nVALUE MISMATCHES:\n")
        for r in bad:
            sys.stderr.write("  %-40s mathilda=%s  wolfram=%s\n"
                             % (r["id"], r["vm"], r["vw"]))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
