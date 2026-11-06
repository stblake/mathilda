#!/usr/bin/env python3
"""
issac97_holdoutc_table.py -- classify the SEALED ISSAC '97 holdout C
(experiments 156-165) against the TRUE answers.

    # run Mathilda on all ten files, save its outputs, classify, write a table
    python3 benchmarks/issac97_holdoutc_table.py --run --mathilda ./Mathilda \\
        --save /tmp/hc --md /tmp/hc/HOLDOUT_C.md
    # classify outputs saved earlier (DIR holds r156.mth ... r165.mth)
    python3 benchmarks/issac97_holdoutc_table.py --from /tmp/hc -v
    # only the stored baseline column
    python3 benchmarks/issac97_holdoutc_table.py --baseline-only
    # a subset
    python3 benchmarks/issac97_holdoutc_table.py --run --only 160,162 --mathilda ./Mathilda

Every case is checked against the TRUE answer listed below (and, with its
provenance, in each experiment's header and in ISSAC97_HOLDOUT_C.md), never
against the other column: a case is

    correct   the CHECK value matches the true answer (one unit of tolerance
              in the last checked digit, for a value that sits on a rounding
              edge; integers and residues must match exactly);
    wrong     a CHECK value that does not match (including a malformed one,
              e.g. a complex number where the answer is real);
    none      no answer: timed out, returned unevaluated, FindRoot failed, a
              SKIP, or no CHECK line at all.

Timing (the BENCH line, milliseconds) is compared only where BOTH columns are
correct.  The baseline column is the stored run in issac97c_baseline/r<NNN>.base
(one other CAS, run once when the holdout was authored).

THE CHECK FORMAT.  sig[v, d] prints {sign * leading d digits, decimal exponent};
a truth entry S(t, d) holds the true value t (to more digits than checked) and
the digit count d of that case; exact parts are Python ints.  A row's truth is
the same nested structure as its CHECK value.
"""

import argparse
import json
import math
import os
import re
import subprocess
import sys
import time
from decimal import Decimal, getcontext

getcontext().prec = 80
HERE = os.path.dirname(os.path.abspath(__file__))
BASELINE_DIR = os.path.join(HERE, "issac97c_baseline")
TIMEOUT_RE = re.compile(r"\(\*\s*bench-timeout:\s*(\d+)\s*\*\)")

EXPS = {
    156: ("156-issac97c-illcond-linalg", "illcond_linalg.m", "c1 ill-conditioned linear algebra"),
    157: ("157-issac97c-hard-integrals", "hard_integrals.m", "c2 hard definite integrals"),
    158: ("158-issac97c-slow-series", "slow_series.m", "c3 slowly convergent series"),
    159: ("159-issac97c-coefficients-residues", "coefficients_residues.m", "c4 coefficients and residues"),
    160: ("160-issac97c-special-roots", "special_roots.m", "c5 zeros of special functions"),
    161: ("161-issac97c-polynomial-systems", "polynomial_systems.m", "c6 polynomial systems"),
    162: ("162-issac97c-cancellation-limits", "cancellation_limits.m", "c7 cancellation limits"),
    163: ("163-issac97c-closed-form-integrals", "closed_form_integrals.m", "c8 closed-form integrals"),
    164: ("164-issac97c-operator-eigenvalues", "operator_eigenvalues.m", "c9 operator eigenvalues"),
    165: ("165-issac97c-ode-singularities", "ode_singularities.m", "c10 ODE singularities"),
}


class S:
    """A sig[v, d] check: true value t (decimal string), d digits."""
    def __init__(self, t, d):
        self.t, self.d = Decimal(t), d

    def __repr__(self):
        return "S(%s, %d)" % (self.t, self.d)


TRUTH = {}


def add(label, spec):
    assert label not in TRUTH, label
    TRUTH[label] = spec


# ---- 156 c1 ill-conditioned linear algebra --------------------------------
add("c1.1 Det exact n=30", [S("2.155511909042743239289243e383", 15), 688963510])
add("c1.1 SingularValueList machine n=8", S("952111739.17589770431", 5))
add("c1.1 SingularValueList wp300 n=60", S("6.5002797341016883363988875221890100680344e123", 20))
add("c1.2 Det exact n=40", [S("4.577979307747334601669105e-95", 15), 680396165, 786603587])
add("c1.2 Eigenvalues machine n=30", S("4.6374958257459897152e-9", 5))
add("c1.2 Eigenvalues wp60 n=120", S("3.6755466792997523573723644127290936748325e-20", 20))
add("c1.3 Inverse exact n=30", [S("9.6076792050570582e16", 15), 378033045])
add("c1.3 SingularValueList machine n=12", S("87522898.09438388228", 5))
add("c1.3 SingularValueList wp250 n=150", S("2.0233509621854076544381811529734466381979e93", 20))
add("c1.4 LinearSolve exact n=60", [2305843009213693950, 2147483646])
add("c1.4 Eigenvalues machine n=20", S("2682440.2336704333605", 6))
add("c1.4 Eigenvalues wp60 n=60", S("2949381270885395289.11123505789", 20))
add("c1.5 SingularValueList machine n=32", S("2000001", 6))
add("c1.5 SingularValueList wp60 n=48", S("2000000000000000000000000000001", 20))
add("c1.6 Det exact n=50", [S("5.493218977963533811877472e-133", 15), 892345966, 335996330])
add("c1.6 Eigenvalues machine n=50", S("0.00050074416716795090139", 7))
add("c1.6 Eigenvalues wp60 n=50", S("5.00493805098963979720392692316e-13", 20))
add("c1.7 Det exact n=14", [S("-3.105041021756156326334828", 15), 351202457, 287719374])
add("c1.7 SingularValueList machine n=8", S("7885.9266359339996083", 8))
add("c1.7 SingularValueList wp60 n=20", S("224185454513517586077.459984416", 20))
add("c1.8 LeastSquares exact n=20", [S("9.500000000005000000005249e10", 15), 934833507, 98001])
add("c1.8 LeastSquares machine n=20", S("9500000.050000525", 6))
add("c1.8 LeastSquares wp80 n=20", S("9500000000000000000000000000000.05", 20))
add("c1.9 Det exact n=30", [S("4.646231590025091403636229e766", 15), 788825102])
add("c1.9 SingularValueList machine n=6", S("350717286.23187149429", 5))
add("c1.9 SingularValueList wp400 n=60", S("3.1446619305684671804581132319806878328733e205", 20))
add("c1.10 LinearSolve exact n=50", S("5868.37134401954916514387168845", 20))
add("c1.10 LinearSolve machine n=60", S("8081580.43332516446719215375233", 5))
add("c1.10 LinearSolve wp60 n=100", S("2261277808846905453744539.11567", 20))

# ---- 157 c2 hard definite integrals ----------------------------------------
for k, t in [(1, "0.164374518416399948150658425158"), (2, "0.624713256427713604289968377817"),
             (4, "998.828034469465151587773873665"), (5, "0.714421505810163771281836190476"),
             (6, "2.45216991381443964446393101662"), (7, "0.536591003574682187526420853985"),
             (9, "2.40411380631918857079947632302"), (10, "1.57079632677179604650584089425")]:
    add("c2.%d NIntegrate machine" % k, S(t, 6))
    add("c2.%d NIntegrate wp40" % k, S(t, 15))
add("c2.3 NIntegrate machine", [S("-1.56898252904292563667809530742", 6), S("2.71755987209473961110709894196", 6)])
add("c2.3 NIntegrate wp40", [S("-1.57079632661351668280789990658", 15), S("2.72069904603716751055632205414", 15)])
add("c2.8 NIntegrate machine", S("6.28287116885660785085983559428e-4", 6))
add("c2.8 NIntegrate wp40", S("6.28318530717644488427169906772e-12", 12))

# ---- 158 c3 slowly convergent series ---------------------------------------
add("c3.1 Sum primes machine", S("0.452247420041065498506543364832", 5))
add("c3.1 MoebiusMu log-zeta wp50", S("0.452247420041065498506543364832248", 25))
add("c3.2 Sum machine", S("0.657973626739290574588966066658", 7))
add("c3.2 Sum Infinity wp30", S("0.657973626739290574588966066658", 20))
for k, t, dm, dw in [(3, "0.582240526465012505902656320160", 8, 20), (4, "7480.97786909167252764630707088", 8, 20),
                     (5, "2.10974280123689197447925719762", 6, 15), (6, "0.924299897222938855959570181360", 8, 20),
                     (8, "4.59987374327233731394301571030", 6, 15), (9, "100.577943338496872490282154286", 8, 20),
                     (10, "0.937548254315843753702574094568", 8, 20)]:
    add("c3.%d NSum machine" % k, S(t, dm))
    add("c3.%d NSum wp40" % k, S(t, dw))
add("c3.7 NSum machine", S("-0.194108935092182049739149244928", 6))
add("c3.7 PolyLog wp40", S("-0.194108935092182049739149244928", 25))

# ---- 159 c4 coefficients and residues --------------------------------------
add("c4.1 SeriesCoefficient exact", [S("1.75155936072122455560e20", 15), 30911638])
add("c4.1 contour machine", S("1.75155936072122455560e20", 6))
add("c4.1 contour Trapezoidal wp40", S("1.75155936072122455560e20", 15))
add("c4.2 SeriesCoefficient exact", [S("9.105434598127360679191946e760", 15), 453374912])
add("c4.2 contour machine n=150", S("2.976345000721957183735135e111", 6))
add("c4.2 contour Trapezoidal wp40", S("9.105434598127360679191946e760", 15))
add("c4.3 SeriesCoefficient exact", [S("1.622236286607641975684111e1529", 15), 476434592])
add("c4.3 contour machine n=150", S("3.178581463702089539995326e113", 6))
add("c4.3 contour Trapezoidal wp40 n=500", S("1.527606893791749160074493e381", 15))
for k, t in [(4, "3.57186931397414815388235508842e-13"), (5, "1.21798130787418619753688753669")]:
    add("c4.%d Residue exact" % k, S(t, 20))
    add("c4.%d contour machine" % k, S(t, 7))
    add("c4.%d contour Trapezoidal wp40" % k, S(t, 20))
add("c4.6 SeriesCoefficient exact", [S("1.8365675328e10", 15), 365675202])
add("c4.6 contour machine", S("1.8365675328e10", 6))
add("c4.6 contour Trapezoidal wp40", S("1.8365675328e10", 11))
add("c4.7 SeriesCoefficient exact", [S("7.995406438144408587511378e626", 15), 27361761])
add("c4.7 contour machine n=60", S("3.477147825267152718277728e86", 6))
add("c4.7 contour Trapezoidal wp40", S("7.995406438144408587511378e626", 15))
add("c4.8 Residue exact", [S("1.208583746537199523332248e-8", 15), 13739, 785096048])
add("c4.8 contour machine", S("1.208583746537199523332248e-8", 7))
add("c4.8 contour Trapezoidal wp40", S("1.208583746537199523332248e-8", 20))
add("c4.9 SeriesCoefficient exact", [S("4.231075653086085499515517e575", 15), 536729466])
add("c4.9 contour machine n=100", S("2.405334743833347895362243e82", 6))
add("c4.9 contour Trapezoidal wp40", S("4.231075653086085499515517e575", 15))
add("c4.10 Residue exact", [S("2.186684242975762143335026e17", 15), 228502357, 371842544])
add("c4.10 contour machine", S("2.186684242975762143335026e17", 6))
add("c4.10 contour Trapezoidal wp40", S("2.186684242975762143335026e17", 20))

# ---- 160 c5 zeros of special functions -------------------------------------
for k, t in [(1, "-2.54264135777352642429380615666"), (3, "4.33323782040642167053239927093"),
             (5, "6.02094890469759665490251152161"), (6, "-2.61072086844414465000153771572"),
             (7, "2.84891782079513985274970473873"), (8, "0.826114765984970336177373238601"),
             (10, "-2.71726282920457410157058066168")]:
    add("c5.%d FindRoot machine" % k, S(t, 8))
    add("c5.%d FindRoot wp40" % k, S(t, 20))
add("c5.2 FindRoot machine", [S("-2.40158510486800288417413977468", 8), S("10.7762995161150708984971033464", 8)])
add("c5.2 FindRoot wp40", [S("-2.40158510486800288417413977468", 20), S("10.7762995161150708984971033464", 20)])
add("c5.4 NSolve machine", S("0.00620338140950646343366027049836", 5))
add("c5.4 NSolve wp60", S("0.00620338140950646343366027049836", 15))
add("c5.9 NSolve machine", [S("-12.0759406304122110376023780407", 5), S("0.718482571514124542165114922967", 5)])
add("c5.9 NSolve wp40", [S("-12.0759406304122110376023780407", 20), S("0.718482571514124542165114922967", 20)])

# ---- 161 c6 polynomial systems ---------------------------------------------
C6 = {1: (3, 2, "0.333333333333333333333333333333", 2, "0.666666666666666666666666666667"),
      2: (4, 4, "2.25", 4, "1.65840319545635272232263624243"),
      3: (4, 4, "-0.531615234375", 4, "0.984693406500657076632830920783"),
      4: (3, 9, "-130", 9, "1.51598022769282058968467151318"),
      5: (2, 16, "3103", 16, "1.79759338472934391991989819789"),
      6: (4, 36, "-0.0001575841326647290296589232", 36, "1.17206247758678450467230451977"),
      7: (4, 14, "0.1375", 21, "1.68372096585234179107896172109"),
      8: (2, 7, "5", 7, "3.23552961063904388209241601064"),
      9: (4, 5, "-357.65625", 5, "6.27469440017220052664355112728"),
      10: (4, 9, "-0.07032260594940866061513273", 9, "3.93580934808760382950002583463")}
for k, (ln, dg, ratio, nsol, mx) in C6.items():
    add("c6.%d GroebnerBasis exact" % k, [ln, dg, S(ratio, 15)])
    add("c6.%d NSolve machine" % k, [nsol, S(mx, 8)])
    add("c6.%d NSolve wp40" % k, [nsol, S(mx, 20)])

# ---- 162 c7 cancellation limits --------------------------------------------
C7 = {1: ("0.0078125", "0.00786110073873594954270274898838", 6, "0.007812500000000000000048828125"),
      2: ("0.072815845483676724860586375874901", "0.0728109999598714357921509401701", 5,
          "0.0728158454836767248605379240589"),
      3: ("-0.5", "-0.604034620028079957818435419433", 6, "-0.507238241365054197127518815315"),
      4: ("0.125", "0.123758547060251570021092889914", 6, "0.12499999999999999999875"),
      5: ("0.051386301650969986602749494316666", "0.051386244424042204412835357969", 6,
          "0.0513863016509699866027494370897"),
      6: ("-0.083333333333333333333333333333333", "-0.0833325055670544577553412563859", 6,
          "-0.0833333333333333333333333333333"),
      7: ("1", "1.00000178871303624146689042945", 7, "1.00000000000000000000000000242"),
      8: ("0.98905599532797255539539565150063", "0.988149496998975588713389726403", 6,
          "0.98905599532797255538632086074"),
      9: ("1", "1.0000016045056037476186910802", 7, "1.00000000000000000000000000235"),
      10: ("-1", "-1.00000715775999526695053169905", 7, "-1.0000000000000000002327585093")}
for k, (lim, p1, dm, p2) in C7.items():
    add("c7.%d Limit exact" % k, S(lim, 20))
    add("c7.%d N machine" % k, S(p1, dm))
    add("c7.%d N wp30" % k, S(p2, 20))

# ---- 163 c8 closed-form integrals ------------------------------------------
for k, t in enumerate(["-1.01494160640965362502120255427", "0.803719849629681710151993020173",
                       "-0.329236162849817068243549448583", "0.0405307333976636291098351317972",
                       "0.523598775598298873077107230547", "-0.260442806300988445540093868790",
                       "0.210714789568552108342911874627", "-0.541161616855569095758001848271",
                       "-2.68811867848189624972422298605", "-0.968946146259369380483634845847"], 1):
    add("c8.%d Integrate exact" % k, S(t, 20))
    add("c8.%d NIntegrate machine" % k, S(t, 7))
    add("c8.%d NIntegrate wp40" % k, S(t, 20))

# ---- 164 c9 operator eigenvalues -------------------------------------------
for lab, t, d in [
        ("c9.1 Eigenvalues machine n=6", "1.18019377358048382524722046390", 8),
        ("c9.1 Eigenvalues wp300 n=60", "1.00019973481797696610152115435", 20),
        ("c9.2 shooting machine", "8.42500692949919857451071877294", 6),
        ("c9.2 shooting wp40", "8.42500692949919857451071877294", 15),
        ("c9.3 Nystrom machine n=20", "0.453572864122606855682199755895", 8),
        ("c9.3 Nystrom wp50 n=50", "0.453572864122606855682199755895", 20),
        ("c9.4 oscillator basis machine N=40", "1.15626707198811325616020504412", 8),
        ("c9.4 oscillator basis wp60 N=100", "1.15626707198811325616020504412", 15),
        ("c9.5 Eigenvalues machine", "0.0379599847847906607432036502309", 8),
        ("c9.5 Eigenvalues wp40", "0.0379599847847906607432036502309", 20),
        ("c9.6 shooting machine", "1.01879297164747108901732478340", 6),
        ("c9.6 shooting wp40", "1.01879297164747108901732478340", 15),
        ("c9.7 Eigenvalues machine N=20", "1.35303016474573485330920174962", 8),
        ("c9.7 Eigenvalues wp60 N=40", "1.35303016474573485330920174962", 20),
        ("c9.8 Eigenvalues machine n=12", "1.12491323229113313151398617253", 8),
        ("c9.8 Eigenvalues wp200 n=40", "1.18367763240231769750996446175", 20),
        ("c9.9 Eigenvalues machine n=8", "0.109050773266525765920701065576", 8),
        ("c9.9 Eigenvalues wp200 n=30", "0.0103729919730295071056672534457", 20),
        ("c9.10 shooting machine", "9.75332212475071491068952390873", 6),
        ("c9.10 shooting wp40", "9.75332212475071491068952390873", 15)]:
    add(lab, S(t, d))

# ---- 165 c10 ODE singularities ---------------------------------------------
for lab, t, d in [
        ("c10.1 NDSolve pole machine", "2.61557120988237384514980859223", 6),
        ("c10.1 NDSolve pole wp40", "2.61557120988237384514980859223", 15),
        ("c10.2 NDSolve pole machine", "1.26277886619324886642763675001", 6),
        ("c10.2 NDSolve pole wp40", "1.26277886619324886642763675001", 15),
        ("c10.3 NDSolve pole machine", "2.84089482764194948767064961441", 6),
        ("c10.3 NDSolve pole wp40", "2.84089482764194948767064961441", 15),
        ("c10.4 NDSolve stiff machine", "40.9452029284009547381334925019", 5),
        ("c10.4 NDSolve stiff wp32", "40.9452029284009547381334925019", 12),
        ("c10.5 NDSolve BVP machine eps=1/100", "1.33597477915866113864683884112", 6),
        ("c10.5 NDSolve BVP wp32 eps=1/1000", "1.33359306918159219890529130807", 12),
        ("c10.6 NDSolve chaotic machine T=50", "-7.06153226337260248572175221451", 5),
        ("c10.6 NDSolve chaotic wp40 T=100", "9.65823343087484059218820188529", 15),
        ("c10.7 NDSolve chaotic machine T=20", "-0.725120577945468811698443108149", 5),
        ("c10.7 NDSolve chaotic wp40 T=40", "0.595491880642409614465877717234", 15),
        ("c10.8 NDSolve singular start machine", "77.4157788604427077180591019652", 6),
        ("c10.8 NDSolve singular start wp40", "77.4157788604427077180591019652", 15),
        ("c10.9 NDSolve machine", "1.33910862429232758493332543647", 6),
        ("c10.9 NDSolve wp40", "1.33910862429232758493332543647", 15),
        ("c10.10 NDSolve backward machine L=8", "0.367061551548078427747792113176", 6),
        ("c10.10 NDSolve backward wp40 L=20", "0.367061551548078427747792113176", 15)]:
    add(lab, S(t, d))


def exp_of(label):
    """'c7.3 Limit exact' -> 162."""
    fam = int(label[1:label.index(".")])
    return 155 + fam


# ---- parsing --------------------------------------------------------------
def parse_output(text):
    """{label: {"ms": float or None, "check": str or None, "skip": bool}}"""
    rows = {}
    for line in text.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3 or parts[0] not in ("BENCH", "CHECK", "SKIP"):
            continue
        r = rows.setdefault(parts[1], {"ms": None, "check": None, "skip": False})
        if parts[0] == "BENCH":
            try:
                r["ms"] = float(parts[2])
            except ValueError:
                pass
        elif parts[0] == "CHECK":
            r["check"] = parts[2].strip()
        else:
            r["skip"] = True
    return rows


def parse_value(s):
    """InputForm of a nested list of integers -> Python lists; a string -> str; else None."""
    s = s.strip()
    if s.startswith('"'):
        return s.strip('"')
    if not re.fullmatch(r"[\s{},\-0-9]+", s):
        return None
    try:
        return json.loads(s.replace("{", "[").replace("}", "]"))
    except ValueError:
        return None


def match(spec, val):
    if isinstance(spec, S):
        if not (isinstance(val, list) and len(val) == 2 and all(isinstance(x, int) for x in val)):
            return False
        m, e = val
        got = Decimal(m).scaleb(e - spec.d + 1)
        t = spec.t
        if t == 0:
            return m == 0
        ulp = Decimal(1).scaleb(int(math.floor(abs(t).log10())) - spec.d + 1)
        return abs(got - t) <= ulp
    if isinstance(spec, int):
        return isinstance(val, int) and val == spec
    if isinstance(spec, list):
        return isinstance(val, list) and len(val) == len(spec) and all(match(a, b) for a, b in zip(spec, val))
    raise TypeError(spec)


def classify(row, spec):
    """-> (status, ms): status in correct / wrong / none."""
    if row is None or row["skip"] or row["check"] is None:
        return "none", None
    v = parse_value(row["check"])
    if isinstance(v, str):
        return "none", None
    if v is None:
        return "wrong", row["ms"]
    return ("correct" if match(spec, v) else "wrong"), row["ms"]


# ---- running Mathilda -----------------------------------------------------
def file_timeout(path):
    try:
        with open(path, encoding="utf-8") as fh:
            m = TIMEOUT_RE.search(fh.read(4096))
        return int(m.group(1)) if m else 1800
    except OSError:
        return 1800


def run_mathilda(binary, num, save):
    d, mfile, _ = EXPS[num]
    cwd = os.path.join(HERE, d)
    tmo = file_timeout(os.path.join(cwd, mfile))
    t0 = time.time()
    try:
        p = subprocess.run([binary, "-file", mfile], cwd=cwd, capture_output=True, text=True,
                           timeout=tmo, stdin=subprocess.DEVNULL)
        out = p.stdout
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        out += "\n# FILE TIMEOUT after %d s\n" % tmo
    print("ran %d in %.0f s" % (num, time.time() - t0), file=sys.stderr)
    if save:
        os.makedirs(save, exist_ok=True)
        with open(os.path.join(save, "r%d.mth" % num), "w") as fh:
            fh.write(out)
    return out


def read(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except OSError:
        return None


# ---- report ---------------------------------------------------------------
def fmt_ms(ms):
    if ms is None:
        return "-"
    return "%.3g" % ms if ms < 100 else "%.0f" % ms


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", action="store_true", help="run Mathilda on the experiments")
    ap.add_argument("--mathilda", default=os.path.join(HERE, "..", "Mathilda"), help="Mathilda binary")
    ap.add_argument("--save", help="directory to save r<NNN>.mth outputs of --run")
    ap.add_argument("--from", dest="from_dir", help="directory holding r<NNN>.mth from an earlier run")
    ap.add_argument("--baseline-only", action="store_true", help="classify only the stored baseline")
    ap.add_argument("--only", help="comma-separated experiment numbers (156-165)")
    ap.add_argument("--md", help="write a markdown table here")
    ap.add_argument("-v", "--verbose", action="store_true", help="list every case")
    a = ap.parse_args()

    nums = sorted(EXPS) if not a.only else [int(x) for x in a.only.split(",")]
    base, mth = {}, {}
    for n in nums:
        t = read(os.path.join(BASELINE_DIR, "r%d.base" % n))
        base.update(parse_output(t) if t else {})
        if a.baseline_only:
            continue
        if a.run:
            t = run_mathilda(os.path.abspath(a.mathilda), n, a.save)
        elif a.from_dir:
            t = read(os.path.join(a.from_dir, "r%d.mth" % n))
        else:
            t = None
        mth.update(parse_output(t) if t else {})
    show_m = not a.baseline_only and (a.run or a.from_dir)

    lines, fam = [], {}
    hdr = "| case | baseline | ms | Mathilda | ms | faster |" if show_m else "| case | baseline | ms |"
    lines += [hdr, "|" + "---|" * (6 if show_m else 3)]
    for label, spec in TRUTH.items():
        n = exp_of(label)
        if n not in nums:
            continue
        bs, bms = classify(base.get(label), spec)
        f = fam.setdefault(n, {"n": 0, "b": {"correct": 0, "wrong": 0, "none": 0},
                               "m": {"correct": 0, "wrong": 0, "none": 0}, "both": 0, "mfast": 0})
        f["n"] += 1
        f["b"][bs] += 1
        row = "| %s | %s | %s |" % (label, bs, fmt_ms(bms))
        if show_m:
            ms_, mms = classify(mth.get(label), spec)
            f["m"][ms_] += 1
            who = ""
            if bs == "correct" and ms_ == "correct" and bms is not None and mms is not None:
                f["both"] += 1
                if mms < bms:
                    f["mfast"] += 1
                who = "Mathilda %.2gx" % (bms / mms) if mms < bms else "baseline %.2gx" % (mms / max(bms, 1e-3))
            row += " %s | %s | %s |" % (ms_, fmt_ms(mms), who)
        lines.append(row)
        if a.verbose:
            print(row)
    # unknown labels (a case the truth table does not know is a harness error)
    for src, rows in (("baseline", base), ("Mathilda", mth)):
        for lab in rows:
            if lab not in TRUTH:
                print("warning: %s case %r has no truth entry" % (src, lab), file=sys.stderr)

    summ = ["", "| family | cases | baseline correct / wrong / none |" +
            (" Mathilda correct / wrong / none | both correct | Mathilda faster |" if show_m else ""),
            "|---|---|---|" + ("---|---|---|" if show_m else "")]
    tot = {"n": 0, "b": [0, 0, 0], "m": [0, 0, 0], "both": 0, "mfast": 0}
    for n in nums:
        f = fam.get(n)
        if not f:
            continue
        b = [f["b"][k] for k in ("correct", "wrong", "none")]
        m = [f["m"][k] for k in ("correct", "wrong", "none")]
        tot["n"] += f["n"]; tot["both"] += f["both"]; tot["mfast"] += f["mfast"]
        tot["b"] = [x + y for x, y in zip(tot["b"], b)]; tot["m"] = [x + y for x, y in zip(tot["m"], m)]
        r = "| %d %s | %d | %d / %d / %d |" % (n, EXPS[n][2], f["n"], *b)
        if show_m:
            r += " %d / %d / %d | %d | %d |" % (*m, f["both"], f["mfast"])
        summ.append(r)
    r = "| **total** | %d | %d / %d / %d |" % (tot["n"], *tot["b"])
    if show_m:
        r += " %d / %d / %d | %d | %d |" % (*tot["m"], tot["both"], tot["mfast"])
    summ.append(r)
    print("\n".join(summ))
    if a.md:
        with open(a.md, "w") as fh:
            fh.write("# ISSAC '97 holdout C (156-165): classification against the true answers\n\n")
            fh.write("\n".join(summ[1:]) + "\n\n" + "\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
