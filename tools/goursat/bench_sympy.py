#!/usr/bin/env python3.11
"""Time the sympy reference (goursat.py) on every GeneralisedGoursat corpus.

    python3.11 tools/goursat/bench_sympy.py [timeout=120] [corpus ...]

Rebuilds the cases exactly as export_corpora.py does (same order, same labels),
runs each in its own process under a wall-clock timeout, and writes
tests/data/generalised_goursat/sympy_<corpus>.tsv with one row per case:

    label  status  seconds  defect  unevaluated  residual

`seconds` covers the reduction only (goursat_integrate / integrate_algebraic), not
the numerical verification, so it is comparable with the Mathilda runner's
AbsoluteTiming of GeneralisedGoursat`Reduce / IntegrateAlgebraic.
"""
import os
import sys
import time
import multiprocessing as mp

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import export_corpora as EC            # noqa: E402  (puts RESEARCH on sys.path)
import sympy as sp                      # noqa: E402
from sympy.parsing.mathematica import parse_mathematica   # noqa: E402
from goursat import goursat_integrate, integrate_algebraic, verify, verify_expr   # noqa: E402

t = EC.t


def corpora():
    def wlc(fname, fam):
        out = []
        for c in EC.wl_list(fname, fam):
            e = parse_mathematica(c["IntegrandWL"]).subs(sp.Symbol('x'), t)
            out.append(dict(c, Integrand=e))
        return out
    return {
        "unit": EC.unit_cases,
        "examples": EC.run_example_cases,
        "literature": EC.literature_cases,
        "welz": EC.welz,
        "stress_2026": lambda: EC.stress_cases(2026),
        "stress_7": lambda: EC.stress_cases(7),
        "stress_11": lambda: EC.stress_cases(11),
        "pseudoelliptic_sqrt": lambda: wlc("goursat_examples.wl", "pseudoelliptic-sqrt"),
        "pseudoelliptic_cbrt": lambda: wlc("goursat_examples_cbrt.wl", "pseudoelliptic-cbrt"),
    }


def work(c, qu):
    try:
        t0 = time.perf_counter()
        if c["Kind"] == "FRp":
            r = goursat_integrate(c["F"], c["R"], t, c["p"])
            secs = time.perf_counter() - t0
            d = verify(r, c["F"], t)
        else:
            r = integrate_algebraic(c["Integrand"], t)
            secs = time.perf_counter() - t0
            d = verify_expr(r, c["Integrand"], t)
        qu.put((r.status, secs, float(d), r.antiderivative.has(sp.Integral),
                str(r.residual).replace("\t", " ")[:200]))
    except Exception as e:      # noqa
        qu.put(("ERROR", time.perf_counter() - t0, None, False, f"{type(e).__name__}: {str(e)[:150]}"))


def run(name, cases, timeout):
    res, running, queue = {}, {}, list(enumerate(cases))
    nproc = max(1, mp.cpu_count() - 2)
    while queue or running:
        while queue and len(running) < nproc:
            i, c = queue.pop(0)
            qu = mp.Queue()
            pr = mp.Process(target=work, args=(c, qu))
            pr.start()
            running[i] = (pr, qu, time.time())
        time.sleep(0.02)
        for i in list(running):
            pr, qu, st = running[i]
            if not qu.empty():
                res[i] = qu.get()
                pr.join()
                del running[i]
            elif time.time() - st > timeout:
                pr.terminate()
                res[i] = ("TIMEOUT", float(timeout), None, False, "")
                del running[i]
            elif not pr.is_alive() and qu.empty():
                res[i] = ("ERROR", 0.0, None, False, "worker died")
                del running[i]
    path = os.path.join(EC.OUT, f"sympy_{name}.tsv")
    with open(path, "w") as f:
        f.write("label\tstatus\tseconds\tdefect\tunevaluated\tresidual\n")
        for i, c in enumerate(cases):
            s, secs, d, un, resid = res[i]
            f.write(f"{c['Label']}\t{s}\t{secs:.4f}\t{'' if d is None else '%.1e' % d}\t{un}\t{resid}\n")
    print(f"{path}: {len(cases)} cases")


if __name__ == "__main__":
    timeout = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    want = sys.argv[2:]
    for name, mk in corpora().items():
        if not want or name in want:
            run(name, mk(), timeout)
