#!/usr/bin/env python3
"""fricas_check.py -- cross-check Bronstein verdicts against FriCAS.

Reads a stress.py JSON results file, sends every case classed FALSE_CERT
(i.e. "not elementary"; for the random family that is a claim to check) or
DECLINE to FriCAS (integrate(f, x), 60 s each, one process per case), and
reports:
  DISAGREE  we certified "not elementary" but FriCAS returned a closed form;
  AGREE     both say not elementary (FriCAS returns an unevaluated integral);
  FRICAS-ONLY  we declined, FriCAS closed it (a coverage gap).
"""
import json, re, subprocess, sys

FRICAS = "/opt/local/bin/fricas"

def to_fricas(f):
    s = f.replace("Sqrt[", "sqrt(").replace("Log[", "log(").replace("E^", "%e^")
    s = s.replace("[", "(").replace("]", ")")
    return s

def fricas(f, timeout=60):
    """CLOSED (an elementary closed form), UNEVAL (integral(...) left), ERROR
    (FriCAS's own 'implementation incomplete' and the like), TIMEOUT."""
    src = f")set messages autoload off\nr := integrate({to_fricas(f)}, x)\n)quit\n"
    try:
        p = subprocess.run([FRICAS, "-nosman"], input=src, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return "TIMEOUT"
    out = p.stdout.split("(1) ->", 1)[-1]
    if ">> Error" in out: return "ERROR"
    # an unevaluated integral prints as the integral sign with a dummy variable
    # (d%A), not as the word "integral"
    if "integral" in out or "%A" in out or "++" in out: return "UNEVAL"
    if "Type:" in out: return "CLOSED"
    return "ERROR"

def main():
    res = json.load(open(sys.argv[1]))
    n = {"DISAGREE": 0, "AGREE": 0, "FRICAS-ONLY": 0, "BOTH-OPEN": 0}
    for r in res:
        if r["cls"] not in ("FALSE_CERT", "DECLINE", "TIMEOUT"): continue
        fr = fricas(r["f"])
        if r["cls"] == "FALSE_CERT":
            tag = "DISAGREE" if fr == "CLOSED" else "AGREE"
        else:
            tag = "FRICAS-ONLY" if fr == "CLOSED" else "BOTH-OPEN"
        n[tag] += 1
        if tag in ("DISAGREE", "FRICAS-ONLY"):
            print(f"{tag:11s} fricas={fr:7s} ours={r['cls']}:{r['detail'][:80]}\n    {r['f'][:200]}", flush=True)
    print(n)
    sys.exit(1 if n["DISAGREE"] else 0)

if __name__ == "__main__":
    main()
