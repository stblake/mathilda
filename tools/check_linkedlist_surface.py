#!/usr/bin/env python3
"""
check_linkedlist_surface.py -- does every head treat a LinkedList as it treats
the List of the same elements, without converting it to one?

A LinkedList (src/linkedlist.h) is a mutable deque of typed or boxed chunks.
The contract (docs/spec/builtins/data-structures.md, LINKED_LIST_DESIGN.md
section 3.2) is that wherever a List works, a LinkedList works too:

  * the same answer: H[ll], with every LinkedList in the result replaced by
    its elements, equals H[Normal[ll]] -- held, so an unevaluated H[ll] is
    never "fixed" by evaluating H[{...}] afterwards. The reference is
    Normal[ll] itself, not the literal it was built from: a typed list's
    Normal is packed, and where the List heads answer a packed and a plain
    List differently (x^0.5 over mixed signs promotes a packed buffer to
    Complex throughout) the LinkedList follows its own Normal;
  * computed ON THE CHUNKS: no List, packed List or NDArray of the list's
    elements is built (Internal`LinkedListConversionCount[] does not move;
    only Normal may materialise a list);
  * the argument untouched (only the mutators change a list).

Nothing here names the heads to check. Every live System` symbol is called in
a handful of shapes over a handful of fixtures, once with plain Lists and once
with the same data as LinkedLists; a call whose List form does not evaluate is
not a list operation and is skipped. A curated EXTRA list adds the multi-
argument constructs the generic shapes cannot spell (iterators, Fold, ...).

Findings, by kind:

    CRASH     the LinkedList call killed the process
    UNEVAL    the List call evaluates, the LinkedList call stays unevaluated
    DIFF      both evaluate, the answers differ
    MSG       a message on one side only
    CONV      the right answer, but the list was materialised to get it
    MUTATED   a non-mutating call changed its LinkedList argument
    LISTSIDE  a DIFF or UNEVAL that is not LinkedList's: the List head itself answers a
              packed List (what Normal[ll] gives) and the same plain List
              differently, and the LinkedList follows one of them. Reported
              for the List head's owner; never gated.

The gate RATCHETS (as check-compile-coverage does): BASELINE_FILE lists the
known findings, one line per (head, kind); a finding not in it fails the run, and a baseline entry no
longer found is reported so it can be deleted (--update rewrites the file).
EXEMPT names heads that are out of scope, each with a reason.

Usage:
    python3 tools/check_linkedlist_surface.py [--mathilda ./Mathilda]
        [--only H1,H2] [-j N] [--update] [--report out.tsv]
"""

import argparse
import concurrent.futures as cf
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BASELINE_FILE = os.path.join(HERE, "data", "linkedlist_surface_baseline.tsv")

# Side effects, interaction, or meaningless on data: never called.
DENY_RE = re.compile(
    r"^(Quit|Exit|Pause|Abort|Interrupt|Throw|Return|Break|Continue|Goto|Label|"
    r"Catch|Sow|Reap|Run|RunProcess|Put|PutAppend|Save|DumpSave|Get|Needs|"
    r"Install|Uninstall|Import|Export|Print|Echo|EchoFunction|EchoTiming|"
    r"Message|Monitor|PrintTemporary|Speak|Beep|Dialog|Input|InputString|"
    r"SystemOpen|CopyToClipboard|Share|Set|SetDelayed|Unset|Clear|ClearAll|"
    r"Remove|Protect|Unprotect|SetAttributes|ClearAttributes|SetOptions|"
    r"ApplyTo|TagSet|TagSetDelayed|TagUnset|UpSet|UpSetDelayed|"
    r"AppendTo|PrependTo|Push|Pop|PopFirst|AddTo|SubtractFrom|TimesBy|DivideBy|"
    r"Increment|Decrement|PreIncrement|PreDecrement|TimeConstrained|"
    r"MemoryConstrained|Compile|LinkedList|Normal|ToExpression|Symbol|Names|"
    r"Definition|Information|SlashCommand|SeedRandom|RandomSeed|Begin|End|"
    r"BeginPackage|EndPackage|SetDirectory|ResetDirectory|SetEnvironment|"
    r"Hold|HoldForm|HoldComplete|HoldPattern|Unevaluated|Evaluate|Trace|"
    r"Timing|AbsoluteTiming|RepeatedTiming|ParallelEvaluate|"
    r"DistributeDefinitions|Locked|FrontEnd.*|Notebook.*|Dynamic.*)$"
    r"|File|Director|Stream|Open|Close|Read|Write|Create|Copy|Rename|Delete.*File"
    r"|Plot|Graphics|Show|Image|Sound|Audio|Video|URL|Socket|Process|Link"
)

# Heads that are out of scope, or whose LinkedList behaviour deliberately
# differs from the List behaviour. One line of reason each.
EXEMPT = {
    "Head": "Head[ll] is LinkedList by design",
    "ListQ": "a LinkedList is not a List (ListQ is False by design)",
    "AtomQ": "a LinkedList is an atom by design",
    "ToString": "prints LinkedList[<n>]: a mutable object has no literal",
    "Hash": "hashes the object identity, by design",
    "ByteCount": "counts the chunk payload, not a List's nodes",
    "SameQ": "identity comparison, by design",
    "UnsameQ": "identity comparison, by design",
    "Order": "identity order, by design",
    "OrderedQ": "OrderedQ[ll] (one atom) is True; contents order via Sort",
    "OwnValues": "fixture symbol bookkeeping",
    "Nothing": "a symbol, not a function of the fixture",
    "MeanFilter": "image filter: numeric List data is a convenience (v0.6740); a LinkedList is not image data",
    "MedianFilter": "image filter: numeric List data is a convenience (v0.6740); a LinkedList is not image data",
    "GaussianFilter": "image filter: numeric List data is a convenience (v0.6740); a LinkedList is not image data",
    "Dilation": "image filter: numeric List data is a convenience (v0.6740); a LinkedList is not image data",
    "Erosion": "image filter: numeric List data is a convenience (v0.6740); a LinkedList is not image data",
    "Opening": "image filter: numeric List data is a convenience (v0.6740); a LinkedList is not image data",
    "Closing": "image filter: numeric List data is a convenience (v0.6740); a LinkedList is not image data",
    "DerivativeFilter": "image filter: numeric List data is a convenience (v0.6740); a LinkedList is not image data",
    "GradientFilter": "image filter: numeric List data is a convenience (v0.6740); a LinkedList is not image data",
    "MorphologicalComponents": "image filter: numeric List data is a convenience (v0.6740); a LinkedList is not image data",
    "Unevaluated": "holds its argument",
    "Equal": "identity of two distinct objects is undecided, by design",
    "Unequal": "identity of two distinct objects is undecided, by design",
    "FullForm": "prints the object",
    "InputForm": "prints the object",
    "TeXForm": "prints the object",
    "Short": "prints the object",
    "Shallow": "prints the object",
    "MatchQ": "{__} patterns do not match a LinkedList (head LinkedList)",
    "RandomChoice": "nondeterministic: two calls never agree",
    "RandomSample": "nondeterministic: two calls never agree",
    "Unique": "makes a fresh symbol every call",
    "Minors": "matrix-only linear algebra: out of scope (user decision)",
    "Orthogonalize": "matrix-only linear algebra (a basis of rows): out of scope (user decision)",
    "PrincipalComponents": "matrix-only linear algebra: out of scope (user decision)",
}

# Matrix-only linear algebra, graphs, ML and plotting are out of scope
# (plan decision 4); they are listed by prefix rather than one by one.
EXEMPT_RE = re.compile(
    r"Graph|Vertex|Edge|Tree|Cluster|Classif|Predict|Learn|Distribution|"
    r"Matrix(?!Q)|Eigen|Inverse$|Det$|LinearSolve|Decomposition|"
    r"SingularValue|NullSpace|RowReduce|MatrixRank|Kronecker|Chart|Histogram"
)

# Fixtures: (name, Mathilda expression for the List form).
FIXTURES = [
    ("int", "{3, 1, 2, 2}"),
    ("real", "{1.5, -2., 3., 0.5}"),
    ("boxed", "{a, b, a, c}"),
    ("str", "{\"b\", \"a\", \"c\", \"a\"}"),
]

SHAPES = ["h[x]", "h[x, 2]", "h[ff, x]", "h[x, ff]", "h[x, EvenQ]",
          "h[x, x]", "h[{x, x}]"]

# Constructs the shapes above cannot spell. `x` is the fixture; `n` is a
# LinkedList (or List) of rows, `nn` one of LinkedLists (or Lists).
EXTRA = [
    "Table[i^2, {i, x}]", "Table[{i, j}, {i, x}, {j, 2}]", "Sum[i, {i, x}]",
    "Product[i, {i, x}]", "Module[{s = 0}, Do[s += i, {i, x}]; s]",
    "Array[ff, Length[x]]", "Fold[ff, 0, x]", "FoldList[ff, 0, x]",
    "MapThread[ff, {x, x}]", "Inner[Times, x, x, Plus]", "Outer[ff, x, x]",
    "Pick[x, {True, False, True, False}]", "Position[x, 2]",
    "Cases[x, _Integer]", "Cases[x, _Real]", "DeleteCases[x, 2]",
    "x /. 2 -> 5", "x /. a -> z", "Replace[x, 2 -> 5, {1}]", "f @@ x",
    "Apply[f, {x}, {1}]", "f @@@ nn", "MapIndexed[ff, x]", "MapAt[ff, x, 2]",
    "Module[{s = {}}, Scan[(s = Append[s, #]) &, x]; s]", "Level[x, {1}]",
    "Partition[x, 2]", "Partition[x, 2, 1]", "Riffle[x, 0]",
    "Thread[ff[x, 1]]", "Thread[x -> 1]", "AllTrue[x, EvenQ]",
    "AnyTrue[x, EvenQ]", "NoneTrue[x, EvenQ]", "SelectFirst[x, EvenQ]",
    "LengthWhile[x, OddQ]", "TakeWhile[x, OddQ]", "FirstPosition[x, 2]",
    "SortBy[x, Minus]", "GatherBy[x, EvenQ]", "SplitBy[x, EvenQ]",
    "Dimensions[n]", "Transpose[n]", "Flatten[n]", "Flatten[nn]",
    "Total[n]", "Total[n, 2]", "Length /@ nn", "Map[Total, nn]",
    "Depth[nn]", "Dimensions[nn]", "MatrixQ[n]", "Count[nn, _Integer, 2]",
    "Position[nn, 2]", "Max[nn]", "Join[x, x]", "Join[n, n, 2]",
    "Dot[x, x]", "Dot[n, n]", "Norm[x]", "x . x", "StringJoin[x]",
    "StringRiffle[x, \",\"]", "Association[Thread[x -> 1]]",
    "AssociationThread[x, x]", "Counts[x]", "FixedPoint[Most, x, 2]",
    "Nest[Rest, x, 2]", "NestList[Rest, x, 2]", "Select[x, EvenQ, 1]",
    "Map[ff, x, {1}]", "Map[ff, nn, {2}]", "ReplacePart[x, 2 -> z]",
    "Extract[x, {{1}, {2}}]", "Take[x, 2]", "Drop[x, 2]", "Reverse[x]",
    "Mean[x]", "Median[x]", "Variance[x]", "Accumulate[x]",
    "Differences[x]", "Sort[x]", "Union[x]", "Tally[x]", "Max[x]",
    "Min[x]", "MinMax[x]", "Clip[x, {0, 2}]", "Rescale[x]",
]

PRELUDE = r"""
(* The List side is Normal of the same LinkedList: the representation (packed
   or not) a LinkedList's elements have, so a packed/unpacked difference of the
   List heads themselves is never charged to LinkedList. *)
mk[l_] := Which[TrueQ[$LL], LinkedList[l], TrueQ[$UNPACKED], Developer`FromPackedArray[Normal[LinkedList[l]]],
  True, Normal[LinkedList[l]]];
fix[] := (ClearAll[ff];
  x = mk[$X];
  n = mk[{{1, 2}, {3, 4}, {5, 6}}];
  nn = mk[{mk[{1, 2}], mk[{3, x0}], mk[{}]}]);
SetAttributes[dn, HoldAll];
(* Association keys are atomic to ReplaceAll (assoc_struct.h), so an
   association is first opened into AssocN[rules] -- on both sides alike --
   and a LinkedList used as a key (Counts[{x, x}]) is then normalised too. *)
dn[r_] := Hold[r] //. {a_Association :> With[{v = AssocN @@ Normal[a]}, v /; True],
   l_LinkedList :> With[{v = Normal[l]}, v /; True]};
shape[r_] := r /. l_LinkedList :> With[{v = LLk[Length[l]]}, v /; True];
R[id_, e_] := Module[{c0, c1, r, m, orig, mut, hd},
  fix[]; orig = If[TrueQ[$LL], Normal[x], x];
  c0 = Internal`LinkedListConversionCount[];
  m = False;
  r = TimeConstrained[Quiet[Check[e, m = True; $Msg]], 4, $TimedOut];
  If[m, r = TimeConstrained[Quiet[e], 4, $TimedOut]];
  c1 = Internal`LinkedListConversionCount[];
  mut = If[TrueQ[$LL], Normal[x] =!= orig, False];
  hd = ToString[Head[r]];
  WriteString["stdout", id, "\t", hd, "\t", If[m, "1", "0"], "\t", c1 - c0,
    "\t", If[mut, "1", "0"], "\t", If[FreeQ[r, _LinkedList], "0", "1"], "\t",
    StringTake[ToString[dn[Evaluate[r]], InputForm], UpTo[600]], "\n"]];
SetAttributes[R, HoldRest];
"""


def run(cmd, timeout):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return p.stdout
    except subprocess.TimeoutExpired as e:
        out = e.stdout
        return out.decode() if isinstance(out, bytes) else (out or "")


def live_heads(mathilda):
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as f:
        f.write('Scan[WriteString["stdout", #, "\\n"] &, Names["System`*"]]\n')
        path = f.name
    names = run([mathilda, "-file", path], 120).split()
    os.unlink(path)
    return sorted(h for h in set(names)
                  if re.fullmatch(r"[A-Z][A-Za-z0-9]*", h) and not DENY_RE.search(h))


def calls(heads, extra=True):
    """[(id, head, fixture_name, fixture_expr, call_text)]"""
    out = []
    for h in heads:
        for k, s in enumerate(SHAPES):
            for fn, fx in FIXTURES:
                out.append((f"{h}#{k}#{fn}", h, fn, fx, s.replace("h[", h + "[", 1)))
    if extra:
        for k, c in enumerate(EXTRA):
            h = re.match(r"[A-Za-z]+", c).group(0)
            for fn, fx in FIXTURES[:3]:
                out.append((f"{h}#x{k}#{fn}", h, fn, fx, c))
    return out


def program(cs, ll, unpacked=False):
    lines = [f"$LL = {'True' if ll else 'False'}; $UNPACKED = {'True' if unpacked else 'False'};", PRELUDE]
    for cid, _h, _fn, fx, c in cs:
        lines.append(f'$X = {fx}; R["{cid}", {c}];')
    return "\n".join(lines) + "\n"


def parse(text):
    res = {}
    for line in text.splitlines():
        p = line.split("\t", 6)
        if len(p) == 7 and "#" in p[0]:
            res[p[0]] = {"head": p[1], "msg": p[2] == "1", "conv": int(p[3] or 0),
                         "mut": p[4] == "1", "ll": p[5] == "1", "val": p[6]}
    return res


def run_batch(mathilda, cs, ll, unpacked=False):
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as f:
        f.write(program(cs, ll, unpacked))
        path = f.name
    got = parse(run([mathilda, "-file", path], 60 + 6 * len(cs)))
    os.unlink(path)
    missing = [c for c in cs if c[0] not in got]
    if missing and len(cs) > 1:          # a crash ends the batch: one by one
        for c in missing:
            got.update(run_batch(mathilda, [c], ll, unpacked) or {})
    elif missing:
        got[missing[0][0]] = {"head": "CRASH", "msg": False, "conv": 0,
                              "mut": False, "ll": False, "val": ""}
    return got


def run_all(mathilda, cs, ll, jobs, chunk=120, unpacked=False):
    res = {}
    parts = [cs[i:i + chunk] for i in range(0, len(cs), chunk)]
    with cf.ThreadPoolExecutor(jobs) as ex:
        for got in ex.map(lambda p: run_batch(mathilda, p, ll, unpacked), parts):
            res.update(got)
    return res


def evaluated(r, head):
    return r and r["head"] not in (head, "CRASH", "$TimedOut", "TimeConstrained")


def classify(cs, plain, linked, unpacked=None):
    rows = []
    for cid, h, fn, _fx, call in cs:
        if h in EXEMPT or EXEMPT_RE.search(h):
            continue
        p, l = plain.get(cid), linked.get(cid)
        if not evaluated(p, h):
            continue                      # not a list operation on this fixture
        if not l:
            continue
        if l["head"] == "CRASH":
            kind = "CRASH"
        elif (l["head"] == h and l["ll"]) or l["head"] == "$TimedOut":
            u = (unpacked or {}).get(cid)
            kind = "LISTSIDE" if u and u["val"] != p["val"] else "UNEVAL"
        elif l["val"] != p["val"]:
            u = (unpacked or {}).get(cid)
            kind = "LISTSIDE" if u and u["val"] != p["val"] else "DIFF"
        elif l["msg"] != p["msg"]:
            kind = "MSG"
        elif l["mut"]:
            kind = "MUTATED"
        elif l["conv"]:
            kind = "CONV"
        else:
            continue
        rows.append((kind, h, cid, call, fn, p["val"], l["val"]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mathilda", default="./Mathilda")
    ap.add_argument("--only", default="")
    ap.add_argument("-j", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--update", action="store_true",
                    help="rewrite the baseline to the current findings")
    ap.add_argument("--report", default="",
                    help="write every finding (with values) to this TSV")
    a = ap.parse_args()

    heads = a.only.split(",") if a.only else live_heads(a.mathilda)
    cs = calls(heads, extra=not a.only)
    if a.only:
        cs += [c for c in calls([], extra=True) if c[1] in heads]
    print(f"check-linkedlist: {len(heads)} heads, {len(cs)} calls", file=sys.stderr)
    plain = run_all(a.mathilda, cs, False, a.j)
    todo = [c for c in cs if evaluated(plain.get(c[0]), c[1])]
    linked = run_all(a.mathilda, todo, True, a.j)
    rows = classify(cs, plain, linked)
    # A DIFF whose List side itself depends on packing is the List head's bug.
    diff_ids = {r[2] for r in rows if r[0] in ("DIFF", "UNEVAL")}
    diffs = [c for c in cs if c[0] in diff_ids]
    if diffs:
        unpacked = run_all(a.mathilda, diffs, False, a.j, unpacked=True)
        rows = classify(cs, plain, linked, unpacked)

    if a.report:
        with open(a.report, "w") as f:
            f.write("kind\thead\tid\tcall\tfixture\tlist\tlinkedlist\n")
            for r in sorted(rows):
                f.write("\t".join(r) + "\n")

    listside = sorted({r[1] for r in rows if r[0] == "LISTSIDE"})
    if listside:
        print("check-linkedlist: List heads that answer packed and plain Lists differently "
              "(not gated): " + ", ".join(listside))
    rows = [r for r in rows if r[0] != "LISTSIDE"]
    found = {f"{r[1]}\t{r[0]}" for r in rows}
    if a.update or not os.path.exists(BASELINE_FILE):
        os.makedirs(os.path.dirname(BASELINE_FILE), exist_ok=True)
        with open(BASELINE_FILE, "w") as f:
            f.write("# head\tkind -- known LinkedList surface gaps "
                    "(tools/check_linkedlist_surface.py --update)\n")
            for k in sorted(found):
                f.write(k + "\n")
        print(f"check-linkedlist: baseline written, {len(found)} findings "
              f"({len({k.split(chr(9))[0] for k in found})} heads)")
        return 0
    with open(BASELINE_FILE) as f:
        base = {l.rstrip("\n") for l in f if l.strip() and not l.startswith("#")}
    if a.only:
        base = {b for b in base if b.split("\t")[0] in heads}
    new, fixed = sorted(found - base), sorted(base - found)
    for k in fixed:
        print("FIXED (delete from baseline):", k.replace("\t", "  "))
    for k in new:
        r = next(r for r in rows if f"{r[1]}\t{r[0]}" == k)
        print(f"NEW {r[0]}: {r[3]}  [{r[4]}]\n    List:       {r[5][:200]}\n"
              f"    LinkedList: {r[6][:200]}")
    print(f"check-linkedlist: {len(found)} findings, {len(new)} new, "
          f"{len(fixed)} fixed, baseline {len(base)}")
    return 1 if new else 0


if __name__ == "__main__":
    sys.exit(main())
