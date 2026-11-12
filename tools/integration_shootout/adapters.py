"""adapters.py -- run a batch of corpus cases through one system (or one
Mathilda Integrate method) and return one raw record per case.

Every adapter makes the system verify ITS OWN answer: it differentiates the
antiderivative, evaluates D[F] - f at the case's three points to >= 30 digits in
its own arithmetic, and reports per point "ok" (|D[F] - f| <= 10^-20 max(1, |f|)),
"bad" or "na" (did not evaluate).  Nothing is translated back between systems.

Raw record: {"id", "status", "time", "nums": [ok|bad|na]*, "cert": bool,
             "text": answer (InputForm-ish, truncated), "detail"}
with status one of ANSWER (the driver classifies it ELEM / SPECIAL / WRONG),
UNEVAL, TIMEOUT, ERROR.

Process model.  Mathilda, wolframscript and sympy start slowly, so they take a
batch of cases per process, each case under the system's own per-case limit
(TimeConstrained / signal.alarm) and the whole process under a deadline; a case
missing from a batch's output (a crash, or a hang the in-system limit did not
catch) is rerun alone, and a lone run that still produces nothing is TIMEOUT
or ERROR.  FriCAS and Maxima start in well under a second and have no reliable
in-system time limit, so they run one process per case under a process timeout.
"""
import json
import os
import re
import shutil
import subprocess
import tempfile
import time

import translate as T

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TOL_EXP = 20
VERIFY = 20          # seconds of in-system verification per point, on top of the timeout
TRUNC = 4000

MATHILDA = os.path.join(REPO, "Mathilda")
WOLFRAM = (shutil.which("wolframscript")
           or "/Applications/Mathematica.app/Contents/MacOS/wolframscript")
FRICAS = shutil.which("fricas") or "/opt/local/bin/fricas"
MAXIMA = shutil.which("maxima") or "/usr/local/bin/maxima"
SYMPY_PY = os.environ.get("SHOOTOUT_SYMPY_PYTHON", "/usr/local/bin/python3.11")


def available(system):
    exe = {"mathilda": MATHILDA, "wolfram": WOLFRAM, "fricas": FRICAS,
           "maxima": MAXIMA, "sympy": SYMPY_PY}[system]
    if not (exe and os.access(exe, os.X_OK)):
        return False
    if system == "sympy":
        return subprocess.run([SYMPY_PY, "-c", "import sympy"], capture_output=True).returncode == 0
    return True


def run_proc(cmd, timeout, stdin=None, cwd=REPO):
    t0 = time.time()
    try:
        p = subprocess.run(cmd, input=stdin, capture_output=True, text=True, timeout=timeout, cwd=cwd)
        return p.stdout, p.stderr, False, time.time() - t0
    except subprocess.TimeoutExpired as e:
        dec = (lambda s: s.decode(errors="replace") if isinstance(s, bytes) else (s or ""))
        return dec(e.stdout), dec(e.stderr), True, time.time() - t0


def parse_rlines(out):
    """@@R<TAB>id<TAB>status<TAB>time<TAB>nums<TAB>cert<TAB>text"""
    recs = {}
    for ln in out.splitlines():
        if not ln.startswith("@@R\t"):
            continue
        p = ln.split("\t", 6) + [""] * 7
        nums = [s for s in p[4].split(",") if s]
        try:
            t = float(p[3])
        except ValueError:
            t = None
        recs[p[1]] = {"id": p[1], "status": p[2], "time": t, "nums": nums,
                      "cert": p[5] == "True", "text": p[6][:TRUNC], "detail": ""}
    return recs


# --------------------------------------------------------------------------
# Wolfram-Language batch (Mathilda and wolframscript share the script)
# --------------------------------------------------------------------------
WL_TMPL = r"""
shootoutRun[id_String, f_, pts_List, tmo_] := Module[{t0, res, cert = False, st, nums = {}, txt, t, r0, garbage},
  t0 = AbsoluteTime[];
  res = TimeConstrained[Check[r0 = INTEGRATE, cert = True; r0, {Integrate::nonelem}], tmo, $Aborted];
  t = AbsoluteTime[] - t0;
  st = Which[res === $Aborted, "TIMEOUT", !FreeQ[res, Integrate], "UNEVAL", True, "ANSWER"];
  (* an "answer" carrying a failure token or a package-private symbol is garbage
     reaching the user (2026-11-12: $Failed + Trager`Private`tragerFail[..]): it
     fails its check outright, so it counts as WRONG, a hard failure *)
  garbage = st === "ANSWER" && (! FreeQ[res, $Failed | $Aborted] ||
    ! FreeQ[res, s_Symbol /; StringContainsQ[Context[s], "Private`"]]);
  If[garbage, nums = {"bad"}];
  If[st === "ANSWER" && ! garbage,
    nums = Table[Module[{d, fv},
        d = TimeConstrained[Quiet[N[(D[res, x] - f) /. x -> p, 40]], VERIFY, $Aborted];
        fv = Quiet[N[f /. x -> p, 40]];
        Which[!NumberQ[d] || !NumberQ[fv], "na",
              Abs[d] <= 10^-TOL Max[1, Abs[fv]], "ok", True, "bad"]], {p, pts}]];
  txt = If[res === $Aborted, "", ToString[res, InputForm]];
  txt = StringReplace[txt, {"\t" -> " ", "\n" -> " "}];
  Print["@@R\t", id, "\t", st, "\t", ToString[N[t]], "\t", StringRiffle[nums, ","], "\t",
    ToString[cert], "\t", StringTake[txt, UpTo[TRUNC]]]];
"""


def wl_script(cases, integrate_call, timeout):
    head = (WL_TMPL.replace("INTEGRATE", integrate_call).replace("VERIFY", str(VERIFY))
            .replace("TOL", str(TOL_EXP)).replace("TRUNC", str(TRUNC)))
    body = "".join(f'shootoutRun["{c["id"]}", {c["f"]}, {c["points"]}, {timeout}];\n' for c in cases)
    return head + body


def mathilda_call(method):
    if method == "Automatic":
        return "Integrate[f, x]"
    return f'Integrate[f, x, Method -> "{method}"]'


def run_wl_batch(exe, cases, call, timeout, tmpdir, startup):
    fd, path = tempfile.mkstemp(suffix=".m", dir=tmpdir)
    with os.fdopen(fd, "w") as fh:
        fh.write(wl_script(cases, call, timeout))
    deadline = startup + len(cases) * (timeout + 3 * VERIFY + 5)
    out, err, to, wall = run_proc([exe, "-file", path], deadline)
    os.unlink(path)
    return parse_rlines(out), err


def batched(runner, cases, timeout, batch):
    """Run in batches; rerun every case missing from its batch alone."""
    recs = {}
    for i in range(0, len(cases), batch):
        chunk = cases[i:i + batch]
        got, _ = runner(chunk)
        recs.update(got)
        for c in chunk:
            if c["id"] in recs:
                continue
            got1, err = runner([c])
            recs[c["id"]] = got1.get(c["id"]) or {
                "id": c["id"], "status": "ERROR" if err.strip() and "TimeConstrained" not in err else "TIMEOUT",
                "time": None, "nums": [], "cert": False, "text": "", "detail": err.strip()[-300:]}
    return recs


def run_mathilda(cases, method, timeout, tmpdir, batch=8):
    call = mathilda_call(method)
    return batched(lambda cs: run_wl_batch(MATHILDA, cs, call, timeout, tmpdir, 20), cases, timeout, batch)


def run_wolfram(cases, timeout, tmpdir, batch=20):
    return batched(lambda cs: run_wl_batch(WOLFRAM, cs, "Integrate[f, x]", timeout, tmpdir, 30),
                   cases, timeout, batch)


# --------------------------------------------------------------------------
# sympy
# --------------------------------------------------------------------------
SYMPY_WORKER = r'''
import json, signal, sys, time
import sympy as sp
from sympy import *
x = sp.Symbol("x")
class TO(Exception): pass
def alarm(*a): raise TO()
signal.signal(signal.SIGALRM, alarm)
TOL = sp.Integer(10) ** -%(tol)d
cases = json.load(open(sys.argv[1])); tmo = int(sys.argv[2])
for c in cases:
    f = eval(c["sympy"])
    st, nums, txt = "ANSWER", [], ""
    t0 = time.time()
    try:
        signal.alarm(tmo)
        r = sp.integrate(f, x)
        signal.alarm(0)
    except TO:
        st = "TIMEOUT"
    except Exception as e:
        signal.alarm(0); st = "ERROR"; txt = type(e).__name__ + ": " + str(e)[:200]
    t = time.time() - t0
    if st == "ANSWER":
        txt = sp.srepr(r) if False else str(r)
        if r.has(sp.Integral):
            st = "UNEVAL"
        else:
            for p in c["points"]:
                try:
                    signal.alarm(%(verify)d)
                    pv = sp.Rational(p)
                    d = sp.N((sp.diff(r, x) - f).subs(x, pv), 40)
                    fv = sp.N(f.subs(x, pv), 40)
                    signal.alarm(0)
                    if not (d.is_number and fv.is_number) or d.has(sp.nan, sp.zoo, sp.oo):
                        nums.append("na")
                    else:
                        nums.append("ok" if abs(d) <= TOL * max(1, abs(fv)) else "bad")
                except Exception:
                    signal.alarm(0); nums.append("na")
    txt = txt.replace("\t", " ").replace("\n", " ")[:%(trunc)d]
    print("@@R\t%%s\t%%s\t%%.3f\t%%s\tFalse\t%%s" %% (c["id"], st, t, ",".join(nums), txt), flush=True)
'''


def run_sympy(cases, timeout, tmpdir, batch=10):
    worker = os.path.join(tmpdir, "sympy_worker.py")
    if not os.path.exists(worker):
        with open(worker, "w") as fh:
            fh.write(SYMPY_WORKER % {"tol": TOL_EXP, "verify": VERIFY, "trunc": TRUNC})

    def runner(cs):
        items = []
        for c in cs:
            items.append({"id": c["id"], "sympy": T.emit(c["tree"], "sympy"), "points": c["points_list"]})
        fd, path = tempfile.mkstemp(suffix=".json", dir=tmpdir)
        with os.fdopen(fd, "w") as fh:
            json.dump(items, fh)
        out, err, to, wall = run_proc([SYMPY_PY, worker, path, str(timeout)],
                                      10 + len(cs) * (timeout + 3 * VERIFY + 5))
        os.unlink(path)
        return parse_rlines(out), err
    return batched(runner, cases, timeout, batch)


# --------------------------------------------------------------------------
# FriCAS (one process per case)
# --------------------------------------------------------------------------
def fricas_float_abs(s):
    """|value| of a FriCAS Float printed by unparse(...::InputForm): float(m, e, 2)."""
    from fractions import Fraction
    s = s.strip()
    m = re.fullmatch(r"float\((-?\d+),\s*(-?\d+),\s*(\d+)\)", s)
    if m:
        return abs(Fraction(int(m.group(1))) * Fraction(int(m.group(3))) ** int(m.group(2)))
    try:
        return abs(Fraction(s))
    except (ValueError, ZeroDivisionError):
        return None


def run_fricas_one(c, timeout, tmpdir):
    from fractions import Fraction
    fd, outfile = tempfile.mkstemp(suffix=".txt", dir=tmpdir)
    os.close(fd)
    os.unlink(outfile)
    f = T.emit(c["tree"], "fricas")
    L = [")set output algebra off", ")set messages type off", ")set messages autoload off",
         f'tf : TextFile := open("{outfile}", "output")',
         f"f := {f}",
         ")set messages time on",
         "r := integrate(f, x)",
         ")set messages time off",
         'writeLine!(tf, concat("RES ", unparse(r::InputForm)))', "flush(tf)"]
    for p in c["points_list"]:
        L.append(f"c := complexNumeric(eval(D(r, x) - f, x, {p}), 40)")
        L.append(f"cf := complexNumeric(eval(f, x, {p}), 40)")
        L.append('writeLine!(tf, concat(["NUM ", unparse(convert(sqrt(norm(c)))@InputForm), " ", '
                 'unparse(convert(sqrt(norm(cf)))@InputForm)]))')
        L.append("flush(tf)")
    L += ["close!(tf)", ")quit"]
    out, err, to, wall = run_proc([FRICAS, "-nosman"], timeout + 3 * VERIFY, stdin="\n".join(L) + "\n")
    lines = open(outfile).read().splitlines() if os.path.exists(outfile) else []
    if os.path.exists(outfile):
        os.unlink(outfile)
    m = re.search(r"Time: .*?= ([\d.]+) sec", out)
    t = float(m.group(1)) if m else round(wall, 3)
    rec = {"id": c["id"], "time": t, "nums": [], "cert": False, "text": "", "detail": ""}
    res = [ln[4:] for ln in lines if ln.startswith("RES ")]
    if not res:
        if to:
            rec.update(status="TIMEOUT", time=None)
        else:
            em = re.search(r">> (?:Error detected within library code:\s*)?(.*?)\n", out)
            rec.update(status="ERROR", detail=(em.group(1).strip() if em else out.strip()[-200:]))
        return rec
    r = res[0]
    rec["text"] = r[:TRUNC]
    if r == "r" or r.startswith("["):
        rec.update(status="ERROR", detail="integrate failed" if r == "r" else "list of results")
        return rec
    if T.uneval_q("fricas", r):
        rec["status"] = "UNEVAL"
        return rec
    rec["status"] = "ANSWER"
    for ln in lines:
        mm = re.match(r"NUM (\S+) (\S+)$", ln)
        if mm:
            dv, fv = fricas_float_abs(mm.group(1)), fricas_float_abs(mm.group(2))
            if dv is None or fv is None:
                rec["nums"].append("na")
            else:
                rec["nums"].append("ok" if dv <= Fraction(1, 10 ** TOL_EXP) * max(1, fv) else "bad")
    while len(rec["nums"]) < len(c["points_list"]):
        rec["nums"].append("na")
    return rec


# --------------------------------------------------------------------------
# Maxima (one process per case)
# --------------------------------------------------------------------------
def run_maxima_one(c, timeout, tmpdir):
    f = T.emit(c["tree"], "maxima")
    pts = "[" + ", ".join(c["points_list"]) + "]"
    script = (
        "display2d:false$ fpprec:40$ ratprint:false$ assume(x>0)$\n"
        f"f: {f}$\n"
        "t0: elapsed_real_time()$\n"
        "r: errcatch(integrate(f, x))$\n"
        "t1: elapsed_real_time()$\n"
        'if r = [] then print("@@ERR") else (r: first(r), print("@@TIME", t1 - t0), '
        'print("@@RES", string(r)), '
        f"for p in {pts} do (d: errcatch(bfloat(cabs(subst(x = p, diff(r, x) - f)))), "
        "fv: errcatch(bfloat(cabs(subst(x = p, f)))), "
        'if d = [] or fv = [] then print("@@NUM na na") else print("@@NUM", first(d), first(fv))))$\n'
    )
    # stdin at EOF: a sign question ("Is x - 1 positive?") must end the run, not wait
    out, err, to, wall = run_proc([MAXIMA, "--very-quiet", "--batch-string=" + script], timeout + 3 * VERIFY, stdin="")
    rec = {"id": c["id"], "time": None, "nums": [], "cert": False, "text": "", "detail": ""}
    if to and not re.search(r"^@@RES ", out, re.M):
        rec["status"] = "TIMEOUT"
        return rec
    m = re.search(r"^@@TIME\s+(\S+)", out, re.M)
    rec["time"] = float(m.group(1)) if m else round(wall, 3)
    m = re.search(r"^@@RES\s+(.*)$", out, re.M)
    if not m or re.search(r"^@@ERR", out, re.M):
        rec.update(status="ERROR", detail=(err or out).strip()[-200:])
        return rec
    r = m.group(1).strip()
    rec["text"] = r[:TRUNC]
    if T.uneval_q("maxima", r):
        rec["status"] = "UNEVAL"
        return rec
    rec["status"] = "ANSWER"
    from fractions import Fraction
    for mm in re.finditer(r"^@@NUM\s+(\S+)\s+(\S+)\s*$", out, re.M):
        a, b = mm.group(1), mm.group(2)
        try:
            dv = Fraction(a.replace("b", "e")) if "b" in a or re.fullmatch(r"-?[\d.]+", a) else None
            fv = Fraction(b.replace("b", "e")) if "b" in b or re.fullmatch(r"-?[\d.]+", b) else None
        except (ValueError, ZeroDivisionError):
            dv = fv = None
        if dv is None or fv is None:
            rec["nums"].append("na")
        else:
            rec["nums"].append("ok" if abs(dv) <= Fraction(1, 10 ** TOL_EXP) * max(1, abs(fv)) else "bad")
    return rec


def run_per_case(fn, cases, timeout, tmpdir):
    return {c["id"]: fn(c, timeout, tmpdir) for c in cases}


def run_system(system, cases, timeout, tmpdir, method=None):
    if system == "mathilda":
        return run_mathilda(cases, method, timeout, tmpdir)
    if system == "wolfram":
        return run_wolfram(cases, timeout, tmpdir)
    if system == "sympy":
        return run_sympy(cases, timeout, tmpdir)
    if system == "fricas":
        return run_per_case(run_fricas_one, cases, timeout, tmpdir)
    if system == "maxima":
        return run_per_case(run_maxima_one, cases, timeout, tmpdir)
    raise ValueError(system)
