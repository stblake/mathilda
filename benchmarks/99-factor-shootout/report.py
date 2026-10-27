#!/usr/bin/env python3
"""Turn results/<tag>-<engine>.json into Markdown tables.

    python3 report.py <tag> [engine ...]      -> Markdown on stdout

Per family: one row per case with the time per engine, Mathilda's ratio to
the fastest correct engine, and Mathilda's status.  A cell shows the time,
or CHECK-FAIL / TIMEOUT / UNSUPPORTED / SKIP / ERROR.  The fastest correct
engine in each row is in bold.

Final correctness is decided here, because it needs every engine: an answer
whose reconstruction passed in run.py is still CHECK-FAIL when its count of
non-constant factors is below the largest count any verified engine reached
(each verified count is a lower bound on the true irreducible count, so the
maximum is the best certificate available; PARI and the reference engine are
always in that set when they ran).  A CHECK-FAIL's time is discarded.
"""
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ORDER = ["mathilda", "reference", "pari", "maxima", "sympy"]
SLOW = 1.2


def fmt(t):
    if t < 0.001:
        return f"{t * 1e6:.0f} µs"
    if t < 1:
        return f"{t * 1e3:.1f} ms"
    return f"{t:.2f} s"


def main():
    tag = sys.argv[1]
    res = {}
    for path in sorted(glob.glob(os.path.join(HERE, "results", f"{tag}-*.json"))):
        res[os.path.basename(path)[len(tag) + 1:-5]] = json.load(open(path))
    engines = sys.argv[2:] or [e for e in ORDER if e in res] + sorted(set(res) - set(ORDER))
    corpus = json.load(open(os.path.join(HERE, "corpus.json")))

    # Final status per (engine, case), after the count comparison.
    final, truth, conflicts = {}, {}, []
    for c in corpus["cases"]:
        recs = {e: res[e].get(c["id"]) for e in engines if res[e].get(c["id"])}
        counts = {e: r["count"] for e, r in recs.items() if r["status"] == "OK"}
        truth[c["id"]] = max(counts.values()) if counts else None
        if len(set(counts.values())) > 1:
            conflicts.append((c["id"], counts))
        for e, r in recs.items():
            st, t = r["status"], r["t"]
            if st == "OK" and r["count"] < truth[c["id"]]:
                st = "CHECK-FAIL"
            if st != "OK":
                t = None
            final[(e, c["id"])] = (st, t, r)

    print(f"# Factor shootout — `{tag}`\n")
    print("Times are the best of 3 in-engine runs after a warm-up (the warm-up alone when it "
          "exceeds 5 s). **Bold** is the fastest correct engine. `ratio` is Mathilda / best "
          "correct engine. CHECK-FAIL = wrong answer (reconstruction failed, or fewer "
          "irreducible factors than another verified engine found); its time is discarded.\n")
    slow, bad, ratios = [], [], {}
    for fam in corpus["families"]:
        print(f"## {fam}\n")
        print("| case | size | " + " | ".join(engines) + " | ratio | Mathilda |")
        print("|---|---|" + "---:|" * len(engines) + "---:|---|")
        for c in [c for c in corpus["cases"] if c["family"] == fam]:
            cells = {e: final.get((e, c["id"]), ("—", None, None)) for e in engines}
            ok = [t for st, t, _ in cells.values() if t is not None]
            best = min(ok) if ok else None
            row = []
            for e in engines:
                st, t, _ = cells[e]
                if t is None:
                    row.append(st)
                else:
                    row.append(f"**{fmt(t)}**" if t == best else fmt(t))
            mst, mt, _ = cells.get("mathilda", ("—", None, None))
            ratio = ""
            if mt is not None and best:
                r = mt / best
                ratios.setdefault(fam, []).append(r)
                ratio = f"{r:.2f}"
                if r > SLOW:
                    slow.append((fam, c["id"], c["size"], r, mt, best))
            elif "mathilda" in engines and mst != "—":
                bad.append((fam, c["id"], c["size"], mst))
            print(f"| {c['id']} | {c['size']} | " + " | ".join(row) + f" | {ratio} | {mst} |")
        print()

    print("## Summary\n")
    print("### Per family (Mathilda / best correct engine)\n")
    print("| family | cases | slower than 1.2× | geo-mean ratio | worst ratio | not OK |")
    print("|---|---:|---:|---:|---:|---:|")
    for fam in corpus["families"]:
        rs = ratios.get(fam, [])
        n = len([c for c in corpus["cases"] if c["family"] == fam])
        nbad = len([b for b in bad if b[0] == fam])
        gm = f"{math.exp(sum(map(math.log, rs)) / len(rs)):.2f}" if rs else "—"
        worst = f"{max(rs):.2f}" if rs else "—"
        print(f"| {fam} | {n} | {len([r for r in rs if r > SLOW])} | {gm} | {worst} | {nbad} |")
    print(f"\n### Cases where Mathilda is slower than the best engine by more than {SLOW}×\n")
    if slow:
        print("| family | case | size | ratio | Mathilda | best |")
        print("|---|---|---|---:|---:|---:|")
        for fam, cid, size, r, mt, best in slow:
            print(f"| {fam} | {cid} | {size} | {r:.2f} | {fmt(mt)} | {fmt(best)} |")
    else:
        print("None.")
    print("\n### Mathilda CHECK-FAIL / TIMEOUT / UNSUPPORTED / ERROR\n")
    if bad:
        print("| family | case | size | status |")
        print("|---|---|---|---|")
        for fam, cid, size, st in bad:
            print(f"| {fam} | {cid} | {size} | {st} |")
    else:
        print("None.")
    if conflicts:
        print("\n### Verified engines that disagree on the factor count\n")
        print("The lower count is under-factored (reported CHECK-FAIL above).\n")
        for cid, counts in conflicts:
            print(f"- {cid}: " + ", ".join(f"{e} {n}" for e, n in counts.items()))


if __name__ == "__main__":
    main()
