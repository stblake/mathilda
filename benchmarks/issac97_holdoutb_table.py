#!/usr/bin/env python3
"""issac97_holdoutb_table.py -- classify ISSAC '97 holdout B (experiments 146-155).

Holdout B is SEALED: 100 problems (ten families of ten) run once, on the final
code, to measure overfitting to the ten challenge problems (126-135) and to
holdout A (136-145).  Problems, true answers and provenance: ISSAC97_HOLDOUT_B.md.

Every case is classified against the TRUE answer, never against the other
column:

    CORRECT    the CHECK value equals the true check value
    WRONG      a value was returned and it is not the true one
    NO ANSWER  "unevaluated" / "timeout" / "failed" / no CHECK line
    SKIP       the head is absent (SKIP line)

and, where both systems are CORRECT, by speed (Mathilda faster or slower than
the baseline column, with the ratio).  The baseline column was recorded once
and is stored in issac97b_baseline/<experiment>.base (the raw tagged output),
so the final run needs only Mathilda.

Usage:
    python3 benchmarks/issac97_holdoutb_table.py --run            # run Mathilda on 146-155, then classify
    python3 benchmarks/issac97_holdoutb_table.py --from DIR       # classify saved Mathilda outputs DIR/<exp>.out
    python3 benchmarks/issac97_holdoutb_table.py --baseline-only  # the stored baseline column against the truth
    options: --mathilda PATH (default ./Mathilda), --only 146,150, --md OUT.md, --save DIR (keep raw outputs)
"""
import argparse, os, re, subprocess, sys, time

BENCH = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BENCH)
BASEDIR = os.path.join(BENCH, "issac97b_baseline")
EXPS = ["146-issac97b-illcond-linalg", "147-issac97b-hard-integrals", "148-issac97b-slow-series",
        "149-issac97b-coefficients-residues", "150-issac97b-special-roots", "151-issac97b-polynomial-systems",
        "152-issac97b-cancellation-limits", "153-issac97b-closed-form-integrals",
        "154-issac97b-operator-eigenvalues", "155-issac97b-ode-singularities"]
FAMILY = {"146": "1 ill-conditioned linear algebra", "147": "2 hard definite integrals", "148": "3 slowly convergent series",
          "149": "4 exact coefficients and residues", "150": "5 zeros of special functions", "151": "6 polynomial systems",
          "152": "7 cancellation limits", "153": "8 closed-form definite integrals", "154": "9 operator eigenvalues",
          "155": "10 ODE singularities"}

# label -> the CHECK value a correct answer prints (InputForm), derived from the
# true answers in ISSAC97_HOLDOUT_B.md with the same sig[v, d] the .m files use.
TRUTH = {
    'b1.1 Eigenvalues machine': '{41552, 9}',
    'b1.1 Eigenvalues wp200': '{10527054388019842647, 69}',
    'b1.2 SingularValueList machine': '{23840, 10}',
    'b1.2 SingularValueList wp150': '{23704575934741660873, 59}',
    'b1.3 LinearSolve exact': '{{121163222332800, 8}, 525390625, 795045632}',
    'b1.3 LinearSolve wp100': '{12116322233280017375, 8}',
    'b1.3 LinearSolve machine': '{151319, 2}',
    'b1.4 SingularValueList machine': '{99569, -6}',
    'b1.4 SingularValueList wp100': '{10267235976092725403, -47}',
    'b1.5 Eigenvalues machine': '{3103, -2}',
    'b1.5 Eigenvalues wp150': '{75090659817199666518, -3}',
    'b1.6 LeastSquares exact': '{{480676444923389, 3}, 250000000, 891207541}',
    'b1.6 LeastSquares wp60': '{48067644492338907213, 3}',
    'b1.6 LeastSquares machine': '{11906073, 1}',
    'b1.7 Det exact': '1',
    'b1.7 Det wp120': 'True',
    'b1.7 Det machine': 'True',
    'b1.8 Eigenvalues machine': '{1778, -8}',
    'b1.8 Eigenvalues wp120': '{12490356420056808796, -29}',
    'b1.9 Eigenvalues machine': '{7428823, -5}',
    'b1.9 Eigenvalues wp60': '{127207212049089, -37}',
    'b1.10 Det exact': '{{285774516383174, -655}, 0, 710434969}',
    'b1.10 Det wp300': '{28577451638317412094, -655}',
    'b1.10 Det machine': '{85416, -21}',
    'b2.1 NIntegrate machine': '{50407, -1}',
    'b2.1 NIntegrate wp30': '{504067061907, -1}',
    'b2.2 NIntegrate machine': '{999760, 2}',
    'b2.2 NIntegrate wp30': '{999760310630, 2}',
    'b2.3 NIntegrate machine': '{314157, 5}',
    'b2.3 NIntegrate wp30': '{314157265359, 5}',
    'b2.4 NIntegrate machine': '{44649, -1}',
    'b2.4 NIntegrate wp30': '{446489755785, -1}',
    'b2.6 NIntegrate machine': '{5996203, -1}',
    'b2.6 NIntegrate wp30': '{599620322995, -1}',
    'b2.8 NIntegrate machine': '{360133, 0}',
    'b2.8 NIntegrate wp30': '{360133023322, 0}',
    'b2.9 NIntegrate machine': '{62513, -3}',
    'b2.9 NIntegrate wp30': '{625129234764, -3}',
    'b2.10 NIntegrate machine': '{11547, 0}',
    'b2.10 NIntegrate wp30': '{115470053838, 0}',
    'b2.5 NIntegrate machine': '{544868, 1083}',
    'b2.5 NIntegrate machine X=25': '{543788, 269}',
    'b2.5 NIntegrate wp30': '{544868449001, 1083}',
    'b2.7 NIntegrate machine': '{1567, 1}',
    'b2.7 NIntegrate wp20': '{156725, 1}',
    'b3.1 NSum machine': '{9242999, -1}',
    'b3.1 NSum wp30': '{924299897223, -1}',
    'b3.2 NSum machine': '{210974, 0}',
    'b3.2 NSum wp30': '{210974280124, 0}',
    'b3.3 NSum machine': '{107080, 0}',
    'b3.3 NSum wp30': '{10707963268, 0}',
    'b3.4 NSum machine': '{937548, -1}',
    'b3.4 NSum wp30': '{937548254316, -1}',
    'b3.5 NSum machine': '{2404114, 0}',
    'b3.5 NSum wp30': '{240411380632, 0}',
    'b3.6 NSum machine': '{6048986, -1}',
    'b3.6 NSum wp30': '{604898643422, -1}',
    'b3.7 NSum machine': '{1598689, -1}',
    'b3.7 NSum wp30': '{159868903742, -1}',
    'b3.8 NSum machine': '{270581, 0}',
    'b3.8 NSum wp30': '{270580808428, 0}',
    'b3.9 NSum machine': '{136843, 0}',
    'b3.9 NSum wp30': '{136843277762, 0}',
    'b3.10 NSum machine': '{50000, -1}',
    'b3.10 NSum wp30': '{500000000000, -1}',
    'b4.1 SeriesCoefficient exact': '{{240614678640326, 31}, 149727991}',
    'b4.1 contour wp60': '{240614678640326, 31}',
    'b4.1 contour machine n=100': '{190569, 8}',
    'b4.2 SeriesCoefficient exact': '{{453899084277566, -79}, 870274057, 0}',
    'b4.2 contour wp160': '{453899084277566, -79}',
    'b4.2 contour machine n=20': '{152234, -4}',
    'b4.3 Residue exact': '{318864502745108, -54}',
    'b4.3 contour wp40': '{318864502745108, -54}',
    'b4.3 contour machine': '{3189, -54}',
    'b4.4 SeriesCoefficient exact': '{{139946140451000, 11}, 946140451}',
    'b4.4 contour wp80': '{13994614045, 11}',
    'b4.4 contour machine n=100': '{293, 2}',
    'b4.5 SeriesCoefficient exact': '{{611327659767719, 472}, 466468457}',
    'b4.5 contour wp100': '{611327659767719, 472}',
    'b4.5 contour machine n=30': '{16974, 12}',
    'b4.6 SeriesCoefficient exact': '{{-303284129702400, 16}, 29760000}',
    'b4.6 contour wp40': '{-303284129702400, 16}',
    'b4.6 contour machine': '{-303284, 16}',
    'b4.7 Residue exact': '{204822323834388, -20}',
    'b4.7 contour wp50': '{204822323834388, -20}',
    'b4.7 contour machine z^3': '{507286, -2}',
    'b4.8 Coefficient exact': '{{700127388815131, 713}, 388014080}',
    'b4.8 contour wp30': '{700127388815131, 713}',
    'b4.8 contour machine n=40': '{210030, 13}',
    'b4.9 SeriesCoefficient exact': '{{160607260103999, 843}, 229457772}',
    'b4.9 contour wp40': '{160607260103999, 843}',
    'b4.9 contour machine n=30': '{846749, 23}',
    'b4.10 Residue exact': '{-200000000000182, 0}',
    'b4.10 contour wp40': '{-200000000000182, 0}',
    'b4.10 contour machine z^10': '{-200199, 0}',
    'b5.1 FindRoot wp60': '{4420915249799639770159036, 1}',
    'b5.1 FindRoot machine n=100': '{1340648734, 1}',
    'b5.2 FindRoot wp40': '{157001048008319383, -3}',
    'b5.2 FindRoot machine': '{1570010480, -3}',
    'b5.3 FindRoot wp40': '{314080729522507863, 3}',
    'b5.3 FindRoot machine': '{3140807295, 3}',
    'b5.5 FindRoot wp40': '{-177009251046479207, 2}',
    'b5.5 FindRoot machine': '{-177009251, 2}',
    'b5.6 FindRoot wp40': '{709575744573535461, 1}',
    'b5.6 FindRoot machine': '{7095757446, 1}',
    'b5.8 FindRoot wp40': '{999937218866126282, -1}',
    'b5.8 FindRoot machine': '{9999372189, -1}',
    'b5.3 BesselJZero wp30': '{31408072952251, 3}',
    'b5.4 FindRoot wp30': '{14194224809460, 3}',
    'b5.4 FindRoot machine': '{1419422481, 3}',
    'b5.4 ZetaZero wp30': '{14194224809460, 3}',
    'b5.5 AiryAiZero wp30': '{-17700925104648, 2}',
    'b5.7 CountRoots exact': '10',
    'b5.7 NSolve wp40': '{10, {208469081014822569, 1}}',
    'b5.7 FindRoot machine': '{2084690810, 1}',
    'b5.9 CountRoots exact': '2',
    'b5.9 NSolve wp80': '{14142135623731, -21}',
    'b5.9 NSolve machine n=8': '{141, -5}',
    'b5.10 FindRoot wp30': '{{14506161632437, 0}, {18809430001533, 0}}',
    'b5.10 FindRoot machine': '{{14506162, 0}, {18809430, 0}}',
    'b6.1 GroebnerBasis exact': '{6, 32, {-131784879900772, -25}}',
    'b6.1 NSolve wp40': '{32, {333333333333333333, -1}}',
    'b6.1 NSolve machine': '{32, {33333333, -1}}',
    'b6.2 GroebnerBasis exact': '{3, 27, {-350540029708690, -1}}',
    'b6.2 NSolve wp40': '{27, {100377595654634978, 0}}',
    'b6.2 NSolve machine': '{27, {10037760, 0}}',
    'b6.3 GroebnerBasis exact': '{4, 14, {-212054783621800, -1}}',
    'b6.3 NSolve wp40': '{14, {150000000000000000, 0}}',
    'b6.3 NSolve machine': '{14, {15000000, 0}}',
    'b6.4 GroebnerBasis exact': '{3, 20, {100000000000000, 0}}',
    'b6.4 NSolve wp40': '{20, {941270710281753576, 0}}',
    'b6.4 NSolve machine': '{20, {94127071, 0}}',
    'b6.7 GroebnerBasis exact': '{2, 36, {152453775851491, -3}}',
    'b6.7 NSolve wp40': '{36, {133135484579998120, 0}}',
    'b6.7 NSolve machine': '{36, {1331, 0}}',
    'b6.8 GroebnerBasis exact': '{4, 4, {-100000000000000, 0}}',
    'b6.8 NSolve wp40': '{24, {175487766624669276, 0}}',
    'b6.8 NSolve machine': '{24, {17548777, 0}}',
    'b6.9 GroebnerBasis exact': '{5, 9, {177777777777778, 2}}',
    'b6.9 NSolve wp40': '{8, {104195973828311374, 1}}',
    'b6.9 NSolve machine': '{8, {10419597, 1}}',
    'b6.10 GroebnerBasis exact': '{6, 10, {205095001667440, -4}}',
    'b6.10 NSolve wp40': '{10, {-672904761331337508, -2}}',
    'b6.10 NSolve machine': '{10, {-67290476, -2}}',
    'b6.5 GroebnerBasis exact': '{1}',
    'b6.5 NSolve wp40': '{0}',
    'b6.5 NSolve machine': '{0}',
    'b6.6 GroebnerBasis eliminate exact': '{1, 4, 3, 9, -1}',
    'b6.6 NSolve wp40': '{8, {916731117205546876, 0}}',
    'b6.6 NSolve machine': '{8, {91673112, 0}}',
    'b7.1 Limit exact': '{270580808427785, -1}',
    'b7.1 approach wp300': '{270580808428, -1}',
    'b7.1 Neville machine': '{2706, -1}',
    'b7.2 Limit exact': '{-738905609893065, 0}',
    'b7.2 Neville wp600': '{-7389056099, 0}',
    'b7.2 Neville machine': '{-74, 0}',
    'b7.3 Limit exact': '{100000000000000, 0}',
    'b7.3 approach wp200': '{100000000000, 0}',
    'b7.3 Neville machine': '{10000, 0}',
    'b7.4 Limit exact': '{728158454836767, -2}',
    'b7.4 approach wp100': '{728158454837, -2}',
    'b7.4 Neville machine': '{72816, -2}',
    'b7.5 Limit exact': '{333333333333333, -1}',
    'b7.5 approach wp150': '{333333333333, -1}',
    'b7.5 Neville machine': '{33333333, -1}',
    'b7.6 Limit exact': '{100000000000000, 0}',
    'b7.6 approach wp150': '{100000000000, 0}',
    'b7.6 sample machine x=1e-8': '{10000, 0}',
    'b7.7 Limit exact': '{-833333333333333, -2}',
    'b7.7 approach wp150': '{-833333333333, -2}',
    'b7.7 Neville machine': '{-8333333, -2}',
    'b7.8 Limit exact': '{-225694444444444, -3}',
    'b7.8 approach wp250': '{-225694444444, -3}',
    'b7.8 Neville machine': '{-226, -3}',
    'b7.9 Limit exact': '{-187500000000000, 0}',
    'b7.9 Neville wp100': '{-187500000000, 0}',
    'b7.9 Neville machine': '{-1875, 0}',
    'b7.10 Limit exact': '{218357912107657, -1}',
    'b7.10 approach wp100': '{218357912108, -1}',
    'b7.10 Neville machine': '{2183579, -1}',
    'b7.3 Series exact': '{100000000000000, 0}',
    'b7.8 Series exact': '{-225694444444444, -3}',
    'b8.1 Integrate exact': '{272198261287950, -1}',
    'b8.1 NIntegrate wp40': '{272198261287950266, -1}',
    'b8.1 NIntegrate machine': '{27219826, -1}',
    'b8.2 Integrate exact': '{-216464646742228, 0}',
    'b8.2 NIntegrate wp40': '{-216464646742227638, 0}',
    'b8.2 NIntegrate machine': '{-216465, 0}',
    'b8.3 Integrate exact': '{183193118835444, 0}',
    'b8.3 NIntegrate wp40': '{183193118835443803, 0}',
    'b8.3 NIntegrate machine': '{18319312, 0}',
    'b8.4 Integrate exact': '{-870057726728316, -1}',
    'b8.4 NIntegrate wp40': '{-870057726728315507, -1}',
    'b8.4 NIntegrate machine': '{-870058, -1}',
    'b8.5 Integrate exact': '{-135290404213892, 0}',
    'b8.5 NIntegrate wp40': '{-135290404213892274, 0}',
    'b8.5 NIntegrate machine': '{-13529040, 0}',
    'b8.6 Integrate exact': '{-548311355616075, -1}',
    'b8.6 NIntegrate wp40': '{-548311355616075479, -1}',
    'b8.6 NIntegrate machine': '{-54831136, -1}',
    'b8.7 Integrate exact': '{-901542677369696, -1}',
    'b8.7 NIntegrate wp40': '{-901542677369695714, -1}',
    'b8.7 NIntegrate machine': '{-901543, -1}',
    'b8.8 Integrate exact': '{386995600539436, -1}',
    'b8.8 NIntegrate wp40': '{386995600539435576, -1}',
    'b8.8 NIntegrate machine': '{38699560, -1}',
    'b8.9 Integrate exact': '{246740110027234, 0}',
    'b8.9 NIntegrate wp40': '{2467401100272340, 0}',
    'b8.9 NIntegrate machine': '{24674011, 0}',
    'b8.10 Integrate exact': '{384373946213433, -1}',
    'b8.10 NIntegrate wp40': '{384373946213432916, -1}',
    'b8.10 NIntegrate machine': '{38437395, -1}',
    'b9.1 secular FindRoot wp40': '{738810809416454997, -1}',
    'b9.1 Nystrom Richardson machine': '{73881081, -1}',
    'b9.2 Gauss-Hermite Nystrom wp30': '{312500000000000, -2}',
    'b9.2 Gauss-Hermite Nystrom machine': '{3125000000, -2}',
    'b9.3 trapezoid Nystrom wp30 N=256': '{7290000000, -1}',
    'b9.3 trapezoid Nystrom machine N=256': '{729000000, -1}',
    'b9.4 Eigenvalues wp50 N=150': '{91798838929666748450, 3}',
    'b9.4 Eigenvalues machine N=500': '{1015240107, 5}',
    'b9.5 Eigenvalues wp60 n=60': '{600000000000000, 1}',
    'b9.5 Eigenvalues machine n=20': '{200000, 1}',
    'b9.6 Fourier matrix wp40': '{-402567795465667873, 1}',
    'b9.6 Fourier matrix machine': '{-402567795466, 1}',
    'b9.6 MathieuCharacteristicA wp30': '{-40256779546567, 1}',
    'b9.6 MathieuCharacteristicA machine': '{-402567795466, 1}',
    'b9.7 harmonic basis wp50 N=100': '{10603620904841828996, 0}',
    'b9.7 harmonic basis machine N=60': '{106036209, 0}',
    'b9.8 Gauss-Legendre Nystrom wp40': '{110210987014802312, -1}',
    'b9.8 Gauss-Legendre Nystrom machine': '{1102109870, -1}',
    'b9.9 Eigenvalues wp40 N=10': '{162028105542010440, -1}',
    'b9.9 Eigenvalues machine N=40': '{11736795, -2}',
    'b9.10 Eigenvalues pencil wp40 N=60': '{100022105340798863, 0}',
    'b9.10 Eigenvalues pencil machine N=400': '{10000051, 0}',
    'b10.1 naive domain end machine': '{2003, 0}',
    'b10.1 reciprocal switch wp40': '{200314735943, 0}',
    'b10.1 reciprocal switch machine': '{200315, 0}',
    'b10.2 EllipticK exact': '{245111152981069, 1}',
    'b10.2 NDSolve quarter period wp40': '{245111152981069, 1}',
    'b10.2 NDSolve quarter period machine': '{245111, 1}',
    'b10.3 NDSolve x(20) wp40': '{13793200, 1}',
    'b10.3 NDSolve x(5) machine': '{-6512, 0}',
    'b10.4 series start + FindRoot wp30': '{6896848619377, 0}',
    'b10.4 series start + FindRoot machine': '{689685, 0}',
    'b10.5 Toepfer wp30': '{3320573362152, -1}',
    'b10.5 Toepfer machine': '{332057, -1}',
    'b10.6 NDSolve y(8) wp50': '{46922076, -8}',
    'b10.6 NDSolve y(3) machine': '{659, -3}',
    'b10.7 naive domain end machine': '{297, 0}',
    'b10.7 energy tail wp30': '{2974477425402, 0}',
    'b10.7 energy tail machine': '{297448, 0}',
    'b10.8 bisection shooting wp30': '{-158807102, 0}',
    'b10.8 bisection shooting machine': '{-158807, 0}',
    'b10.9 naive domain end machine': '{2221, 0}',
    'b10.9 inverse function wp30': '{2221441469079, 0}',
    'b10.9 inverse function machine': '{222144, 0}',
    'b10.10 NDSolve e=0.99 wp30': '{{-905786813, -1}, {140566255, -1}}',
    'b10.10 NDSolve e=0.9 machine': '{{-7147, -1}, {4283, -1}}',
}

NOANSWER = {'"unevaluated"', '"timeout"', '"failed"'}

def norm(s):
    return re.sub(r"\s+", "", s or "")

def parse(text):
    bench, check, skip = {}, {}, {}
    for line in text.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and f[0] == "BENCH":
            try: bench[f[1]] = float(f[2])
            except ValueError: pass
        elif len(f) >= 3 and f[0] == "CHECK":
            check[f[1]] = f[2]
        elif len(f) >= 3 and f[0] == "SKIP":
            skip[f[1]] = f[2]
    return bench, check, skip

def timeout_of(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        m = re.search(r"\(\*\s*bench-timeout:\s*(\d+)\s*\*\)", fh.read(4096))
    return int(m.group(1)) if m else 600

def mfile(exp):
    d = os.path.join(BENCH, exp)
    return d, [f for f in os.listdir(d) if f.endswith(".m")][0]

def run_mathilda(exp, binary, save):
    d, m = mfile(exp)
    t0 = time.time()
    try:
        p = subprocess.run([binary, "-file", m], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           stdin=subprocess.DEVNULL, timeout=timeout_of(os.path.join(d, m)) + 60)
        out = p.stdout.decode("utf-8", "replace")
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or b"").decode("utf-8", "replace")
    print("  %s: %.0f s" % (exp, time.time() - t0), file=sys.stderr)
    if save:
        os.makedirs(save, exist_ok=True)
        with open(os.path.join(save, exp + ".out"), "w") as fh: fh.write(out)
    return out

def klass(label, check, skip):
    if label in skip: return "SKIP"
    v = check.get(label)
    if v is None or norm(v) in NOANSWER: return "NO ANSWER"
    return "CORRECT" if norm(v) == norm(TRUTH[label]) else "WRONG"

def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--run", action="store_true", help="run Mathilda on each experiment")
    g.add_argument("--from", dest="src", help="directory of saved Mathilda outputs <exp>.out")
    g.add_argument("--baseline-only", action="store_true", help="classify the stored baseline column alone")
    ap.add_argument("--mathilda", default=os.environ.get("MATHILDA_BIN", os.path.join(ROOT, "Mathilda")))
    ap.add_argument("--only", help="comma-separated experiment numbers")
    ap.add_argument("--save", help="with --run: keep the raw outputs here")
    ap.add_argument("--md", help="also write the report to this file")
    a = ap.parse_args()
    exps = [e for e in EXPS if not a.only or e[:3] in a.only.split(",")]
    L, rows = [], []
    for exp in exps:
        num = exp[:3]
        with open(os.path.join(BASEDIR, exp + ".base"), encoding="utf-8", errors="replace") as fh:
            bb, bc, bs = parse(fh.read())
        if a.baseline_only:
            mt = ""
        elif a.run:
            mt = run_mathilda(exp, a.mathilda, a.save)
        else:
            with open(os.path.join(a.src, exp + ".out"), encoding="utf-8", errors="replace") as fh: mt = fh.read()
        mb, mcheck, ms = parse(mt)
        labels = [l for l in TRUTH if l.startswith("b%d." % (int(num) - 145))]
        for lab in labels:
            km, kb = klass(lab, mcheck, ms), klass(lab, bc, bs)
            tm = mb.get(lab) if km == "CORRECT" else None
            tb = bb.get(lab) if kb == "CORRECT" else None
            speed = ""
            if tm is not None and tb is not None:
                speed = "faster" if tm < tb else "slower"
            rows.append(dict(num=num, label=lab, km=km, kb=kb, tm=tm, tb=tb, speed=speed,
                             vm=mcheck.get(lab, ""), vb=bc.get(lab, "")))
    L.append("# ISSAC '97 holdout B (experiments 146-155): classification against the true answers\n")
    L.append("| family | cases | Mathilda correct | baseline correct | both correct | Mathilda faster | Mathilda slower | only Mathilda correct | only baseline correct | Mathilda wrong | Mathilda no answer |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    tot = None
    def summary(rs):
        c = lambda f: sum(1 for r in rs if f(r))
        return [len(rs), c(lambda r: r["km"] == "CORRECT"), c(lambda r: r["kb"] == "CORRECT"),
                c(lambda r: r["km"] == r["kb"] == "CORRECT"), c(lambda r: r["speed"] == "faster"),
                c(lambda r: r["speed"] == "slower"), c(lambda r: r["km"] == "CORRECT" and r["kb"] != "CORRECT"),
                c(lambda r: r["kb"] == "CORRECT" and r["km"] != "CORRECT"), c(lambda r: r["km"] == "WRONG"),
                c(lambda r: r["km"] in ("NO ANSWER", "SKIP"))]
    for exp in exps:
        rs = [r for r in rows if r["num"] == exp[:3]]
        L.append("| %s | %s |" % (FAMILY[exp[:3]], " | ".join(str(x) for x in summary(rs))))
    L.append("| **all** | %s |" % " | ".join("**%d**" % x for x in summary(rows)))
    L.append("\n## Cases\n")
    L.append("| case | Mathilda | baseline | Mathilda ms | baseline ms | ratio | speed |")
    L.append("|---|---|---|---|---|---|---|")
    for r in rows:
        ratio = "%.2f" % (r["tm"] / r["tb"]) if r["tm"] is not None and r["tb"] not in (None, 0.0) else ""
        L.append("| %s | %s | %s | %s | %s | %s | %s |" % (r["label"], r["km"], r["kb"],
                 "" if r["tm"] is None else "%.1f" % r["tm"], "" if r["tb"] is None else "%.1f" % r["tb"], ratio, r["speed"]))
    wrong = [r for r in rows if r["km"] == "WRONG"]
    if wrong:
        L.append("\n## Mathilda wrong answers\n")
        L.append("| case | Mathilda check | true check |")
        L.append("|---|---|---|")
        for r in wrong: L.append("| %s | `%s` | `%s` |" % (r["label"], r["vm"], TRUTH[r["label"]]))
    bad = [r for r in rows if r["kb"] != "CORRECT"]
    if bad:
        L.append("\n## Baseline column: cases not answered correctly\n")
        L.append("| case | baseline | baseline check | true check |")
        L.append("|---|---|---|---|")
        for r in bad: L.append("| %s | %s | `%s` | `%s` |" % (r["label"], r["kb"], r["vb"], TRUTH[r["label"]]))
    text = "\n".join(L) + "\n"
    print(text)
    if a.md:
        with open(a.md, "w") as fh: fh.write(text)

if __name__ == "__main__":
    main()
