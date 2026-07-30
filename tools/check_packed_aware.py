#!/usr/bin/env python3
"""
check_packed_aware.py -- does every head with an NDArray fast path opt in to it?

THE FAILURE MODE THIS EXISTS FOR. The packing transparency gate (src/eval.c,
step 2.7) materialises a packed buffer for any head that has not opted in via
src/pack.c's AWARE list. That rule makes a missing opt-in fail SAFE -- the answer
is still right -- and therefore SILENT. Nothing fails, nothing warns, and the
operation quietly runs at one Expr per element.

It has happened four times, at 30x-658x each:

    the 26 linear-algebra heads   LinearSolve of a 90x90 real system
                                  STACK-OVERFLOWED, because the gate materialised
                                  first, linalg_call_has_ndarray then read false,
                                  and a machine-real solve took the exact
                                  fraction-free path
    Nest / NestList / FixedPoint  118x -- Fold was aware, Nest was not
    the structural family         30x-237x (RotateLeft, Join, Partition,
                                  Differences, Riffle, PadLeft, PadRight)
    Fourier                       40x-80x, still open

So: a head that HAS a fast path but is not on the AWARE list is a bug, and this
script is the thing that says so. It reads the source rather than the running
system because that is where dispatch lives.

A head may be absent from AWARE deliberately -- List enforces the no-nesting
invariant, Floor and friends answer with different element heads. Those go in
EXEMPT below WITH A REASON. An unexplained absence is an error.

Usage:  python3 tools/check_packed_aware.py        (exit 1 on a finding)
        make check-packed-aware
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")

# ---------------------------------------------------------------------------
# Heads that are deliberately NOT packed-aware. Every entry needs a reason; the
# reason is the point of the list, so that a future reader can tell "considered
# and rejected" from "never noticed".
# ---------------------------------------------------------------------------
EXEMPT = {
    "List": "enforces the no-nesting invariant -- a buffer must never sit inside "
            "a plain List, which is what keeps the gate an O(argc) top-level scan",
    # Floor / Ceiling / Round / IntegerPart / Sign were exempt here until
    # 2026-07-30 for want of a narrowing float64->int64 kernel category. They
    # have one now (NDUnaryKernel.to_int), are off NOT_AWARE, and answer with the
    # exact Integers the List path does -- so they are no longer exceptions.
    "Im": "projection kernel (to_real), so a real input yields a real 0.0 where "
          "the List gives the exact Integer 0; it needs the narrowing treatment "
          "Floor and friends got, which its to_real category cannot express",
    "Clip": "clamps to the EXACT bounds, so a clipped element comes back Integer",
    "Precision": "Listable here, so the list answers per element and a buffer "
                 "path would answer with one scalar",
    "Accuracy": "as Precision",

    # Surfaced by this script on 2026-07-30 and then REJECTED by a differential
    # sweep: each answers differently for a buffer than for the same value as a
    # plain List. The observed difference is recorded so nobody re-litigates it.
    "Chop": "replaces a negligible element with the EXACT Integer 0, so "
            "Chop[{1.*^-14, 1.}] is {0, 1.} and a float64 buffer makes it "
            "{0., 1.} -- the Floor class again",
    "Quotient": "real in, exact Integer out: Quotient[{1.,2.,3.,4.}, 2] is "
                "{0,1,1,2} and the buffer path gives {0.,1.,1.,2.} "
                "(HPC_IMPROVEMENT_PLAN.md phase 2 narrowing kernels)",
    "LeafCount": "counts a packed list as ONE node (1 against 5 for a "
                 "4-element vector), which would perturb Simplify's complexity "
                 "metric -- the reason it was excluded by design",
    "ConjugateTranspose": "its NDArray path does not handle rank 1 and comes "
                          "back UNEVALUATED as Conjugate[Transpose[v]]; a real "
                          "gap, but fixing it is not a packing change",
    "Series": "has no NDArray handling; the call stays unevaluated either way, "
              "but the packed form prints its argument differently",
    "SortBy": "its NDArray path does not sort -- SortBy[buffer, Abs] returns the "
              "input unchanged. A live gap on the VISIBLE NDArray surface too, "
              "but fixing it is a sort change, not a packing one",
    "ReverseSortBy": "delegates to SortBy, so it inherits that gap",
    "MinimalBy": "flagged only by this script's intra-file call-graph "
                 "propagation from src/sort.c; maximal_minimal_by walks "
                 "data.function.args directly and has no NDArray handling at all",
    "MaximalBy": "as MinimalBy",
}


def read(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return f.read()


def all_sources():
    for dirpath, _dirs, files in os.walk(SRC):
        if os.sep + "external" in dirpath:
            continue
        for fn in files:
            if fn.endswith(".c"):
                yield os.path.join(dirpath, fn)


def builtin_registry():
    """builtin C function name -> registered Mathilda head name."""
    reg = {}
    pat = re.compile(r'symtab_add_builtin\(\s*"([A-Za-z$][A-Za-z0-9$]*)"\s*,\s*'
                     r'([A-Za-z_][A-Za-z0-9_]*)\s*\)', re.S)
    pat_sym = re.compile(r'symtab_add_builtin\(\s*SYM_([A-Za-z0-9_]+)\s*,\s*'
                         r'([A-Za-z_][A-Za-z0-9_]*)\s*\)', re.S)
    for path in all_sources():
        text = read(path)
        for name, fn in pat.findall(text):
            reg[fn] = name
        for name, fn in pat_sym.findall(text):
            reg.setdefault(fn, name)
    return reg


FN_DEF = re.compile(r'^(?:static\s+)?Expr\*\s*([A-Za-z_][A-Za-z0-9_]*)\s*\(', re.M)

# Evidence, inside a builtin's body, that the head has an NDArray fast path.
#
# The first group is the shared dispatch/degrade helpers. The second group is a
# DIRECT test on the tag, and leaving it out was a false negative that cost a
# real one: src/fourier.c tests `data->type == EXPR_NDARRAY` and calls
# machine_path_ndarray, a complete fast path -- but it uses none of the shared
# helpers, so the audit passed while Fourier[packedList] ran 15.7x slower than
# Fourier[NDArray[...]] on identical data. Any way of noticing the tag counts.
DISPATCH_MARKERS = (
    "linalg_call_has_ndarray",
    "ndstruct_call_has_ndarray",
    "ndstruct_delist_repack",
    "ndarray_delist_and_reeval",
    "linalg_delist_and_reeval",
    "is_ndarray(",
    "EXPR_NDARRAY",
    "is_packed_list(",
)


CALL = re.compile(r'\b([A-Za-z_][A-Za-z0-9_]*)\s*\(')


def heads_with_fast_paths(reg):
    """registered head name -> (file, the marker that identified it).

    Attribution is per function, then propagated ONE CALL AT A TIME through the
    file's own static helpers to a fixed point. The propagation is not
    optional: src/fourier.c tests the NDArray tag inside `fourier_compute`, a
    static helper two calls below `builtin_fourier`, so a function-local scan
    credited nothing and Fourier passed the audit while running 15.7x slower on
    a packed list than on the identical visible NDArray.

    Only intra-file calls are followed. A fast path implemented in another
    translation unit is reached through one of the shared helpers in
    DISPATCH_MARKERS, which are matched directly.
    """
    found = {}
    for path in all_sources():
        text = read(path)
        if not any(m in text for m in DISPATCH_MARKERS):
            continue
        bounds = [(m.start(), m.group(1)) for m in FN_DEF.finditer(text)]
        bodies, marked = {}, {}
        for i, (pos, fn) in enumerate(bounds):
            end = bounds[i + 1][0] if i + 1 < len(bounds) else len(text)
            bodies[fn] = text[pos:end]
            for marker in DISPATCH_MARKERS:
                if marker in bodies[fn]:
                    marked[fn] = marker
                    break
        # Propagate: a function that calls a marked function is itself marked.
        changed = True
        while changed:
            changed = False
            for fn, body in bodies.items():
                if fn in marked:
                    continue
                for callee in set(CALL.findall(body)):
                    if callee != fn and callee in marked:
                        marked[fn] = marked[callee]
                        changed = True
                        break
        for fn, marker in marked.items():
            if fn in reg:
                found.setdefault(reg[fn], (os.path.relpath(path, ROOT), marker))
    return found


def kernel_heads():
    """Heads with a registered ND kernel -- implicitly packed_aware."""
    out = set()
    pat = re.compile(r'symtab_set_ndarray_(?:unary|binary)_kernel\(\s*"?([A-Za-z]+)"?')
    # NOT line-anchored: ndkernels_init packs several per line
    # (`REG_U(Floor); REG_U(Ceiling); REG_U(Round);`), and anchoring found only
    # the first of each, which reported the rest as stale NOT_AWARE entries.
    macro = re.compile(r'REG_[UB]\(\s*([A-Za-z0-9_]+)')
    for path in all_sources():
        text = read(path)
        out.update(pat.findall(text))
        out.update(macro.findall(text))
    return out


def strip_comments(text):
    """Remove /* ... */ and // ... so quoted prose in a comment is not mistaken
    for a list entry. It was: a comment inside the NOT_AWARE initialiser reading
    `the "optimisation, not a correctness question" the note below anticipated`
    was parsed as a head named exactly that."""
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


def string_list(text, varname):
    """Parse `static const char* const VARNAME[] = { "A", "B", ... };`."""
    m = re.search(r'\b' + varname + r'\s*\[\s*\]\s*=\s*\{(.*?)\}\s*;', text, re.S)
    if not m:
        return None
    return set(re.findall(r'"([^"]+)"', strip_comments(m.group(1))))


def main():
    pack_c = read(os.path.join(SRC, "pack.c"))
    aware = string_list(pack_c, "AWARE")
    not_aware = string_list(pack_c, "NOT_AWARE")
    int64_ok = string_list(pack_c, "INT64_OK")
    if aware is None or not_aware is None or int64_ok is None:
        print("ERROR: could not parse the AWARE / NOT_AWARE / INT64_OK lists in "
              "src/pack.c -- has their shape changed?", file=sys.stderr)
        return 2

    reg = builtin_registry()
    fast = heads_with_fast_paths(reg)
    kernels = kernel_heads()

    # symtab_set_ndarray_*_kernel sets packed_aware itself, so a head with a
    # registered kernel is aware WITHOUT appearing in the AWARE list -- and is
    # un-aware only by appearing in NOT_AWARE. Both checks below have to reason
    # about the effective set, not the literal list; UnitStep, which is aware
    # purely through its kernel, was reported as a problem when they did not.
    effective_aware = (aware | kernels) - not_aware

    problems = []

    # (1) A head with an NDArray dispatch that is not aware and not exempt.
    for head, (path, marker) in sorted(fast.items()):
        if head in effective_aware or head in EXEMPT:
            continue
        problems.append(
            f"{head}: has an NDArray fast path ({marker} in {path}) but is not in\n"
            f"    src/pack.c's AWARE list, so the gate materialises the buffer "
            f"before it can be used.\n"
            f"    Add it to AWARE, or add it to EXEMPT in this script with a reason.")

    # (2) INT64_OK implies packed_aware; an entry that is not aware does nothing.
    for head in sorted(int64_ok - effective_aware):
        if head in EXEMPT:
            continue
        problems.append(
            f"{head}: is in INT64_OK but not in AWARE. INT64_OK is the narrower "
            f"claim and\n    has no effect without the broader one.")

    # (3) A head cleared by NOT_AWARE with no registered kernel has nothing to
    #     clear. That is harmless (clearing an unset bit is a no-op) and is
    #     sometimes deliberate belt-and-braces, so it is reported as a NOTE
    #     rather than a failure -- but it is worth seeing, because it also
    #     describes an entry that has drifted away from the code it guards.
    notes = [f"{h}: in NOT_AWARE but has no registered ND kernel -- the clear is "
             f"a no-op (defensive, or stale)"
             for h in sorted(not_aware - kernels)]

    print(f"scanned {len(reg)} registered builtins; "
          f"{len(fast)} have an NDArray dispatch, "
          f"{len(kernels)} have an ND kernel")
    print(f"AWARE={len(aware)}  NOT_AWARE={len(not_aware)}  INT64_OK={len(int64_ok)}  "
          f"EXEMPT={len(EXEMPT)}")

    for n in notes:
        print("  note: " + n)

    if problems:
        print("\npacked-aware audit FAILED:\n", file=sys.stderr)
        for p in problems:
            print("  " + p + "\n", file=sys.stderr)
        print(f"{len(problems)} problem(s). See docs/design/packed_arrays.md.",
              file=sys.stderr)
        return 1

    print("OK: every head with an NDArray fast path opts in (or is exempt with a "
          "reason).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
