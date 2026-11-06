#!/usr/bin/env python3
"""ISSAC '97 holdout (experiments 136-145): classify a run against the TRUE answers.

    python3 benchmarks/issac97_holdout_table.py RESULTS_DIR > table.md

RESULTS_DIR holds, per experiment N in 136..145, the raw output of the .m file
in Mathilda as rN.mth and in the baseline system as rN.base (BENCH/CHECK lines,
see benchmarks/harness.m).  Unlike run_all.py, which compares the two systems'
CHECK values with each other, this script compares each system with the
independently established true answer (listed per problem in each .m file's
header and in benchmarks/ISSAC97_HOLDOUT.md), so a CHECK-FAIL says WHICH system
is wrong.

Classes (Mathilda's point of view):
    AHEAD       both correct, Mathilda time / baseline time < 1.5
    SLOWER      both correct, ratio >= 1.5
    INCOMPLETE  Mathilda gave no answer (timeout / unevaluated / failed)
    CHECK-FAIL  a system answered and its check differs from the truth
                (the "who" column names the system(s) that match the truth)
Requires mpmath (for the leading-digit checks of the true values).
"""
import os
import re
import sys
from mpmath import mp, mpf, floor, log10, nint

mp.dps = 60


def sig(v, d):
    """The harness's sig[v, d]: {sign * leading d digits, decimal exponent}."""
    v = mpf(v)
    if v == 0:
        return "{0, 0}"
    a = abs(v)
    e = int(floor(log10(a)))
    m = int(nint(a / mpf(10) ** (e - d + 1)))
    if m == 10 ** d:              # rounding carried into a new digit
        e += 1
        m = int(nint(a / mpf(10) ** (e - d + 1)))
    return "{%d, %d}" % (m if v > 0 else -m, e)


# ---- true values ----------------------------------------------------------
# (label prefix, kind) -> value; kinds are matched against the case label.
T = {}


def num(label, value, d):
    T[label] = sig(value, d)


def raw(label, s):
    T[label] = s


# family 1: machine 4 digits, wp 12
F1 = {
    "h1.1": ("8.763949109211799472e11", "9.191581826323644019411559631e116"),
    "h1.2": ("5.670141114735457244638779e-9", "1.214052974061843949823397917e-42"),
    "h1.3": ("1.147829915846752502509278903e-8", "1.580934514544763948518980013e-63"),
    "h1.4": ("6.9281547456602970221054600808e-7", "4.488013531293049259832894189e-20"),
    "h1.5": ("0.03102806064401001507954023064", "0.007509065981719966651839414854"),
    "h1.6": ("2.0217595347358266054975227987e-6", "1.272072120490893356864979269e-37"),
    "h1.7": ("1.778272859770121495153843015e-8", "6.738206080100924396271878802e-48"),
    "h1.8": ("5.3644058795433677680264085377e-7", "6.158049892052418768707796020e-48"),
    "h1.9": ("1.196723905131748921764892147e6", "1.199177720852744324270695235e32"),
    "h1.10": ("8.6512575126143135188029573691e-10", "5.423169387591532736872088203e-21"),
}
for k, (m, w) in F1.items():
    num(k + " machine", m, 4)
    num(k + " wp", w, 12)

# family 2: machine 6, wp40 15
F2 = {
    "h2.1": ("1.476250472508459122292813e151", "1.390195540965462192941711e2180"),
    "h2.2": ("6.207592136617010347600381e35", "4.880478191980581312448124e698"),
    "h2.3": ("0.5040670619069283719898561",) * 2,
    "h2.4": ("1.291285997062663540407283",) * 2,
    "h2.5": ("1.826649837699537925233829e-9", "2.283395465013097564662795e-869"),
    "h2.6": ("15707963266.94896619231322",) * 2,
    "h2.7": ("0.4464897557846246056092822",) * 2,
    "h2.8": ("999.2042906713881633188507",) * 2,
    "h2.9": ("9.332621544394415268169924e157", "4.023872600770937735437024e2567"),
    "h2.10": ("0.07786967112327907846744718",) * 2,
}
for k, (m, w) in F2.items():
    num(k + " machine", m, 6)
    num(k + " wp", w, 15)

# family 3: machine 6, wp30 14
F3 = {
    "h3.1": "100.5779433384968724902821543", "h3.2": "2.109742801236891974479257198",
    "h3.3": "0.6048986434216303702472659142", "h3.4": "0.1598689037424309717569478703",
    "h3.5": "8.957773692316167741131746530", "h3.6": "1.070796326794896619231321692",
    "h3.7": "-0.1941089350921820497391492449", "h3.8": "38.40676809282178631849374770",
    "h3.9": "0.5222702824645705054936517728", "h3.10": "19.64669717802233418224120349",
}
for k, v in F3.items():
    num(k + " machine", v, 6)
    num(k + " wp", v, 14)

# family 4: exact {sig15, num mod 1e9, den mod 1e9} (pure integer arithmetic);
# contour machine 5 digits (reduced instance), wp 14 (full instance)
F4 = {
    "h4.1": ("{{190945496569613, 952}, 573069977, 1}", "9.168205127076471545e93", "1.9094549656961307771e952"),
    "h4.2": ("{{398492964932431, 1051}, 321023520, 1}", "7.169778242085354729e103", "3.9849296493243108281e1051"),
    "h4.3": ("{{283775735672486, 992}, 276670786, 1}", "6.229683408098898651e97", "2.8377573567248613515e992"),
    "h4.4": ("{{115859142706818, 16}, 270681804, 1}", "3876087157", "119680802813"),
    "h4.5": ("{{-100297809149655, -94}, 111993200, 140464321}", "-4.8231843915270198307517654368e-11",
             "-1.0029780914965471547121571775e-94"),
    "h4.6": ("{{-665628459853324, 449}, 526592864, 1}", "3.9029050826561728283e41", "-6.6562845985332426505e449"),
    "h4.7": ("{{173102400594873, 27}, 633786324, 1}", "8679280", "1.7310240059487250166e27"),
    "h4.8": ("{{238955623730970, 901}, 204849471, 1}", "1.2515277677761546113e89", "2.3895562373096964257e901"),
    "h4.9": ("{{288961131710638, -79}, 297218968, 533203125}", "1.15876444327988522001979320291e-8",
             "2.88961131710637621756341532854e-79"),
    "h4.10": ("{{563484790092564, -2}, 667105165, 104412672}", "0.176197052001953125",
              "0.0563484790092564222472452641429"),
}
for k, (e, m, w) in F4.items():
    raw(k + " exact", e)
    num(k + " contour NIntegrate machine", m, 5)
    num(k + " contour Trapezoidal", w, 14)

# family 5: machine 9, wp50 20
F5 = {
    "h5.1": ("13.406487338144910138498015358651607", "31.050746380090020220242053712643817"),
    "h5.2": ("2.88870192448943012370974812175e-6",) * 2,
    "h5.3": ("2.22786116597166768664029545987e-4",) * 2,
    "h5.4": ("3.0811692773857049463438593532248656e-7",) * 2,
    "h5.5": ("313.37426607752784471969024510190451",) * 2,
    "h5.6": ("108.83616589840977436309799199049782",) * 2,
    "h5.7": ("316.14950153608385672196249252951105", "3167.3374738866212685428555577555482"),
    "h5.8": ("0.00072271758021023132698199952253706034",) * 2,
    "h5.9": ("-38.021008677255254433132468290748645",) * 2,
    "h5.10": ("49.773832477672302181916784678563724",) * 2,
}
for k, (m, w) in F5.items():
    num(k + " FindRoot machine", m, 6 if k in ("h5.2", "h5.3", "h5.4") else 9)
    num(k + " FindRoot wp", w, 20)

# family 6: exact basis structure + ratio (16 digits); NSolve {count, Max Re}
E4 = "{{0, 0, 0, 16}, {0, 0, 1, 15}, {0, 1, 0, 15}, {1, 0, 0, 15}}"
F6 = {
    "h6.1": (4, E4, "0.86965788645256176780912388945606013321", 16, "0.49979737359804185931674311923711859991"),
    "h6.2": (3, "{{0, 0, 27}, {0, 1, 26}, {1, 0, 26}}", "0.01501362839974060569124851539110572881", 27,
             "3.1343815800329427962215604140507932310"),
    "h6.3": (5, "{{0, 0, 0, 0, 32}, {0, 0, 0, 1, 31}, {0, 0, 1, 0, 31}, {0, 1, 0, 0, 31}, {1, 0, 0, 0, 31}}",
             "192.07810134279376134363036323083912860336", 32, "3.9220137669196961228913674877851927363"),
    "h6.4": (12, "{{1, 1, 2, 2}, {1, 1, 2, 2}, {1, 1, 2, 3}, {1, 1, 2, 3}, {1, 1, 2, 3}, {1, 1, 2, 4}, "
                 "{1, 1, 2, 4}, {1, 1, 2, 4}, {1, 1, 2, 5}, {1, 1, 3, 3}, {1, 2, 2, 2}, {2, 1, 2, 2}}",
             None, 16, "0.90169897874995578505044728772243330428"),
    "h6.5": (6, "{{0, 0, 0, 0, 0, 32}, {0, 0, 0, 0, 1, 31}, {0, 0, 0, 1, 0, 31}, {0, 0, 1, 0, 0, 31}, "
                "{0, 1, 0, 0, 0, 31}, {1, 0, 0, 0, 0, 31}}",
             "9.97684590967955719376317908624395601813515443388e-11", 32, "0.31004455812896180309"),
    "h6.6": (11, "{{0, 0, 0, 0, 15}, {0, 0, 0, 2, 12}, {0, 0, 0, 7, 12}, {0, 0, 1, 0, 11}, {0, 0, 1, 6, 12}, "
                 "{0, 0, 3, 6, 13}, {0, 1, 0, 0, 11}, {0, 1, 0, 6, 12}, {0, 1, 2, 6, 12}, {0, 2, 0, 6, 12}, "
                 "{1, 1, 1, 1, 1}}", "-1", 70, "2.1180339887498948482045868343656381177"),
    "h6.7": (2, "{{0, 56}, {1, 55}}", "1.39346228018519819172976340609918087806", 56,
             "1.5993259998540608616306369836527993470"),
    "h6.8": (4, E4, "36.78270028974930448372755064505423192774", 16, "3.3591528937659752676136105651342438089"),
    "h6.9": (4, "{{0, 0, 0, 24}, {0, 0, 1, 23}, {0, 1, 0, 23}, {1, 0, 0, 23}}",
             "7.17251366003849767913656645595062136834", 24, "0.99786418492623171021759669616277129140"),
    "h6.10": (4, "{{0, 0, 0, 36}, {0, 0, 1, 35}, {0, 1, 0, 35}, {1, 0, 0, 35}}",
              "-32.0637409838524774157642215008156258788", 36, "3.3881757126573176537361304837532591831"),
}
for k, (n, ex, ratio, ns, mx) in F6.items():
    raw(k + " GroebnerBasis", "{%d, %s, %s}" % (n, ex, '"none"' if ratio is None else sig(ratio, 16)))
    raw(k + " NSolve machine", "{%d, %s}" % (ns, sig(mx, 8)))
    raw(k + " NSolve wp", "{%d, %s}" % (ns, sig(mx, 18)))

# family 7: exact limit (20 digits), machine point (5; h7.3 4), wp point (20)
F7 = {
    "h7.1": ("1", "1.01673590039224256361398", "1.000000000000000166666666667"),
    "h7.2": (mpf(4096) / 729, "0.7586975708497706115883995", "5.618655692684505410760728417"),
    "h7.3": (mpf(9765625) / 1048576, "0.7999597155712885935357546", "9.313225746154790314690520366"),
    "h7.4": ("-1", "-1.105145831114269734995316", "-1.001000500166708341668055754"),
    "h7.5": (-mp.e ** 2, "-6.150321866214310661305405", "-7.374307488864213536806713516"),
    "h7.6": (mp.e ** (-mpf(1) / 2), "0.6256388323119630929767135", "0.606530659712633423603799535"),
    "h7.7": (mpf(3) / 2, "1.503869110345547681674374", "1.500000000000000000000000375"),
    "h7.8": (mpf(-17) / 2520, "-0.006767975352357744219015611", "-0.006746031746031746031767901235"),
    "h7.9": ("1", "1.051089324812166718686341", "1.000000000000000500000000000000107"),
    "h7.10": (mpf(1) / 64, "2.063705497271428362992129", "0.01562499996571180559310121929"),
}
for k, (l, m, w) in F7.items():
    num(k + " Limit", l, 20)
    num(k + " point machine", m, 4 if k == "h7.3" else 5)
    num(k + " point N", w, 20)

# family 8: exact (20), machine 8, wp30 14
F8 = {
    "h8.1": "-2.1646464674222763830320073930823", "h8.2": "0.23025792994521775240374192579582",
    "h8.3": "0.14251419793571091570872150120965", "h8.4": "-0.96894614625936938048363484584692",
    "h8.5": "0.17282745097458205019574093418642", "h8.6": "0.37931985286416790881897546186484",
    "h8.7": "-6.493939402266829149096022179247", "h8.8": "-0.75128556447474642837483635094466",
    "h8.9": "3.8757845850374775219345393833877", "h8.10": "-0.30765458032881955246584991473604",
}
for k, v in F8.items():
    num(k + " Integrate", v, 20)
    num(k + " NIntegrate machine", v, 8)
    num(k + " NIntegrate wp", v, 14)

# family 9: machine 9, wp40 18 (h9.9: 7 / 14)
F9 = {
    "h9.1": "1.353030164745734853309201749615", "h9.2": "0.5362060346834763502347030920059",
    "h9.3": "1.253314101714016907331673330101", "h9.4": "0.5359573175891977563406853434110",
    "h9.5": "1.424941177795263197428278866833", "h9.6": "1.665886980965023224079412839504",
    "h9.7": "-5.800046020851508419103472514311", "h9.8": "1.060362090484182899647046016693",
    "h9.9": "6.548395306000592515133370778167", "h9.10": "3.509384443552449794992539728544",
}
for k, v in F9.items():
    num(k + " machine", v, 7 if k == "h9.9" else 9)
    num(k + " wp", v, 14 if k == "h9.9" else 18)

# family 10: machine 6 (h10.8 naive 5), wp30 15
F10 = {
    "h10.1": "2.0031473594268847080046109790543", "h10.2": "1.2091995761561452337293855050948",
    "h10.3": "6.8968486193769603754545281871231", "h10.4": "3.6537537362191224246090942804593",
    "h10.5": "0.54935272877527081901868319424682", "h10.6": "0.33205733621519629893718006201058",
    "h10.7": "1.15959526696392836576999205157", "h10.8": "1.9863527074304728134718183899204",
    "h10.10": "-0.24593576445134833519776086248533",
}
for k, v in F10.items():
    num(k + " machine", v, 5 if k == "h10.8" else 6)
    num(k + " wp", v, 15)
r, yr = "1.61713753496866934400101139575", "1.1495671346072722154505435895"
raw("h10.9 inverse function machine", "{%s, %s}" % (sig(r, 6), sig(yr, 6)))
raw("h10.9 inverse function wp", "{%s, %s}" % (sig(r, 15), sig(yr, 15)))


def truth_for(label):
    """Longest key of T that the label starts with (keys use a word prefix:
    'h1.1 machine' matches 'h1.1 eigenvalues machine n=12')."""
    best = None
    head = label.split(" ")[0]
    for k in T:
        kh, rest = k.split(" ", 1)
        if kh != head:
            continue
        words = rest.split(" ")
        if all(w in label.split(" ") or any(t.startswith(w) for t in label.split(" ")) for w in words):
            if best is None or len(k) > len(best):
                best = k
    return T[best] if best else None


def parse(path):
    bench, check = {}, {}
    if not os.path.exists(path):
        return bench, check
    with open(path, errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[0] == "BENCH" and len(p) >= 3:
                bench[p[1]] = float(p[2])
            elif p[0] == "CHECK" and len(p) >= 3:
                check[p[1]] = p[2].strip()
    return bench, check


def main():
    d = sys.argv[1]
    rows, summary = [], {}
    names = {136: "ill-conditioned matrices", 137: "hard definite integrals", 138: "slow series",
             139: "exact coefficients", 140: "special-function roots", 141: "polynomial systems",
             142: "cancellation limits", 143: "polylog integrals", 144: "operator eigenvalues",
             145: "ODE singularities"}
    for n in range(136, 146):
        mb, mc = parse(os.path.join(d, "r%d.mth" % n))
        bb, bc = parse(os.path.join(d, "r%d.base" % n))
        fam = summary.setdefault(n, {"AHEAD": 0, "SLOWER": 0, "INCOMPLETE": 0, "CHECK-FAIL": 0})
        for lab in mc:
            t = truth_for(lab)
            mans, bans = lab in mb, lab in bb
            mok = mans and mc.get(lab) == t
            bok = bans and bc.get(lab) == t
            ratio = ""
            if not mans:
                cls = "INCOMPLETE"
                who = "baseline correct" if bok else ("baseline wrong" if bans else "baseline also none")
            elif not mok or (bans and not bok):
                cls = "CHECK-FAIL"
                who = {(True, False): "Mathilda", (False, True): "baseline",
                       (False, False): "neither"}.get((mok, bok), "")
                if mok and not bans:
                    who = "Mathilda (baseline none)"
                if not mok and not bans:
                    who = "neither (baseline none)"
            elif not bans:
                cls, who = "AHEAD", "baseline none"
            else:
                r = mb[lab] / bb[lab] if bb[lab] > 0 else float("inf")
                ratio = "%.2f" % r
                cls, who = ("SLOWER" if r >= 1.5 else "AHEAD"), ""
            fam[cls] += 1
            mt = "%.1f" % mb[lab] if mans else mc[lab].strip('"')
            bt = "%.1f" % bb[lab] if bans else bc.get(lab, "-").strip('"')
            rows.append((lab, mt, bt, ratio, cls, who, mc.get(lab, ""), bc.get(lab, ""), t))
    print("| case | Mathilda ms | baseline ms | ratio | class | matches truth |")
    print("|---|---|---|---|---|---|")
    for lab, mt, bt, ratio, cls, who, *_ in rows:
        print("| %s | %s | %s | %s | %s | %s |" % (lab, mt, bt, ratio, cls, who))
    print()
    print("| exp | family | AHEAD | SLOWER | INCOMPLETE | CHECK-FAIL |")
    print("|---|---|---|---|---|---|")
    tot = {"AHEAD": 0, "SLOWER": 0, "INCOMPLETE": 0, "CHECK-FAIL": 0}
    for n, s in summary.items():
        print("| %d | %s | %d | %d | %d | %d |" % (n, names[n], s["AHEAD"], s["SLOWER"], s["INCOMPLETE"], s["CHECK-FAIL"]))
        for k in tot:
            tot[k] += s[k]
    print("| | **total** | %d | %d | %d | %d |" % (tot["AHEAD"], tot["SLOWER"], tot["INCOMPLETE"], tot["CHECK-FAIL"]))
    if "-v" in sys.argv:
        print()
        for lab, mt, bt, ratio, cls, who, mcv, bcv, t in rows:
            if cls in ("CHECK-FAIL", "INCOMPLETE"):
                print("%-50s M=%s  B=%s  T=%s" % (lab, mcv, bcv, t))


if __name__ == "__main__":
    main()
