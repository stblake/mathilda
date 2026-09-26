# -*- coding: utf-8 -*-
"""
First-principles stress corpus for Refine[expr, assum].

Each case:
  id       : unique short id
  cat      : category label
  expr     : the full Mathilda expression whose result is under test
             (usually Refine[...]; may be a raw Assuming/Block wrapper)
  expected : a Mathilda expression the result SHOULD be SameQ (===) to,
             derived from mathematical first principles
  baseline : the "unchanged" form; if result === baseline (and != expected)
             the verdict is UNCHANGED (a missed reduction). None to skip.
  note     : the first-principles justification for `expected`
  hyp      : my prior hypothesis  pass | gap | limit | sound  (for the report;
             NOT used by the verdict logic)

Verdict (computed by run_stress.py from the run, never from `hyp`):
  PASS       result === expected
  UNCHANGED  result === baseline and not PASS         (missing feature)
  DIVERGENT  neither                                   (hand-review: wrong vs
                                                        sound-but-incomplete)
  HANG/CRASH process timeout / signal
"""

CASES = [

# =====================================================================
# 1. Radicals & powers   -- sqrt of squares, principal-branch powers
# =====================================================================
{"id":"rad01","cat":"1 Radicals/Powers","expr":"Refine[Sqrt[x^2], x > 0]","expected":"x","baseline":"Sqrt[x^2]","note":"sqrt(x^2)=|x|=x for x>0","hyp":"pass"},
{"id":"rad02","cat":"1 Radicals/Powers","expr":"Refine[Sqrt[x^2], x < 0]","expected":"-x","baseline":"Sqrt[x^2]","note":"sqrt(x^2)=|x|=-x for x<0","hyp":"pass"},
{"id":"rad03","cat":"1 Radicals/Powers","expr":"Refine[Sqrt[x^2], Element[x, Reals]]","expected":"Abs[x]","baseline":"Sqrt[x^2]","note":"sqrt(x^2)=|x| for real x, sign unknown","hyp":"pass"},
{"id":"rad04","cat":"1 Radicals/Powers","expr":"Refine[Sqrt[x^2], x >= 0]","expected":"x","baseline":"Sqrt[x^2]","note":"x>=0 => |x|=x","hyp":"pass"},
{"id":"rad05","cat":"1 Radicals/Powers","expr":"Refine[Sqrt[x^2], x <= 0]","expected":"-x","baseline":"Sqrt[x^2]","note":"x<=0 => |x|=-x (NonPositive)","hyp":"gap"},
{"id":"rad06","cat":"1 Radicals/Powers","expr":"Refine[Sqrt[x^2], x == 0]","expected":"0","baseline":"Sqrt[x^2]","note":"x=0 => sqrt(0)=0","hyp":"pass"},
{"id":"rad07","cat":"1 Radicals/Powers","expr":"Refine[Sqrt[x^2 y^2], x > 0 && y < 0]","expected":"-x y","baseline":"Sqrt[x^2 y^2]","note":"|x||y|=x*(-y)","hyp":"pass"},
{"id":"rad08","cat":"1 Radicals/Powers","expr":"Refine[Sqrt[x^2 y^2], x > 0 && y > 0]","expected":"x y","baseline":"Sqrt[x^2 y^2]","note":"|x||y|=xy","hyp":"pass"},
{"id":"rad09","cat":"1 Radicals/Powers","expr":"Refine[Sqrt[x^2 y^2], Element[x | y, Reals]]","expected":"Abs[x] Abs[y]","baseline":"Sqrt[x^2 y^2]","note":"sqrt(x^2 y^2)=|x||y| for real x,y","hyp":"gap"},
{"id":"rad10","cat":"1 Radicals/Powers","expr":"Refine[(x^2)^r, x > 0]","expected":"x^(2 r)","baseline":"(x^2)^r","note":"(x^2)^r=x^(2r) for x>0","hyp":"pass"},
{"id":"rad11","cat":"1 Radicals/Powers","expr":"Refine[(x^m)^r, x > 0]","expected":"x^(m r)","baseline":"(x^m)^r","note":"(x^m)^r=x^(mr) for x>0","hyp":"pass"},
{"id":"rad12","cat":"1 Radicals/Powers","expr":"Refine[(a^b)^c, -1 < b < 1]","expected":"a^(b c)","baseline":"(a^b)^c","note":"principal branch for -1<b<1","hyp":"pass"},
{"id":"rad13","cat":"1 Radicals/Powers","expr":"Refine[a^p b^p, a > 0 && b > 0]","expected":"(a b)^p","baseline":"a^p b^p","note":"a^p b^p=(ab)^p for a,b>0","hyp":"pass"},
{"id":"rad14","cat":"1 Radicals/Powers","expr":"Refine[(x^3)^(1/3), x >= 0]","expected":"x","baseline":"(x^3)^(1/3)","note":"cube-root of cube for x>=0","hyp":"pass"},
{"id":"rad15","cat":"1 Radicals/Powers","expr":"Refine[Sqrt[-x^2], Element[x, Reals]]","expected":"I Abs[x]","baseline":"Sqrt[-x^2]","note":"sqrt(-x^2)=i|x| for real x","hyp":"pass"},
{"id":"rad16","cat":"1 Radicals/Powers","expr":"Refine[Sqrt[1/x], x > 0]","expected":"x^(-1/2)","baseline":"Sqrt[1/x]","note":"sqrt(1/x)=x^(-1/2) for x>0","hyp":"pass"},
{"id":"rad17","cat":"1 Radicals/Powers","expr":"Refine[(x^4)^(1/2), x > 0]","expected":"x^2","baseline":"(x^4)^(1/2)","note":"(x^4)^(1/2)=x^2 for x>0","hyp":"pass"},
{"id":"rad18","cat":"1 Radicals/Powers","expr":"Refine[Sqrt[x^2], x > 0 && Element[x, Reals]]","expected":"x","baseline":"Sqrt[x^2]","note":"redundant real fact + positivity","hyp":"pass"},

# =====================================================================
# 2. Abs / Sign / Conjugate / Arg
# =====================================================================
{"id":"abs01","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Abs[x], x > 0]","expected":"x","baseline":"Abs[x]","note":"|x|=x for x>0","hyp":"pass"},
{"id":"abs02","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Abs[x], x < 0]","expected":"-x","baseline":"Abs[x]","note":"|x|=-x for x<0","hyp":"pass"},
{"id":"abs03","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Abs[x], Element[x, Reals]]","expected":"Abs[x]","baseline":"Abs[x]","note":"real, sign unknown -> stays (sound)","hyp":"pass"},
{"id":"abs04","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Sign[x], x > 0]","expected":"1","baseline":"Sign[x]","note":"sign(x)=1 for x>0","hyp":"pass"},
{"id":"abs05","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Sign[x], x < 0]","expected":"-1","baseline":"Sign[x]","note":"sign(x)=-1 for x<0","hyp":"pass"},
{"id":"abs06","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Sign[x], x >= 0]","expected":"Sign[x]","baseline":"Sign[x]","note":"x>=0: sign is 1 for x>0, 0 at 0 -> not constant, stays (sound)","hyp":"pass"},
{"id":"abs07","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Abs[x], x <= 0]","expected":"-x","baseline":"Abs[x]","note":"|x|=-x for x<=0 (NonPositive)","hyp":"gap"},
{"id":"abs08","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Abs[x y], x > 0 && y > 0]","expected":"x y","baseline":"Abs[x y]","note":"|xy|=xy for x,y>0","hyp":"pass"},
{"id":"abs09","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Sign[x y], x > 0 && y < 0]","expected":"-1","baseline":"Sign[x y]","note":"xy<0 => sign=-1","hyp":"gap"},
{"id":"abs10","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Arg[x], x > 0]","expected":"0","baseline":"Arg[x]","note":"arg of positive real =0","hyp":"gap"},
{"id":"abs11","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Arg[x], x < 0]","expected":"Pi","baseline":"Arg[x]","note":"arg of negative real =Pi","hyp":"gap"},
{"id":"abs12","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Conjugate[x], Element[x, Reals]]","expected":"x","baseline":"Conjugate[x]","note":"conj of real = itself","hyp":"pass"},
{"id":"abs13","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Conjugate[a + b I], Element[a | b, Reals]]","expected":"a - I b","baseline":"Conjugate[a + b I]","note":"conj(a+bi)=a-bi for real a,b","hyp":"pass"},
{"id":"abs14","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Abs[x^2 + 1], Element[x, Reals]]","expected":"1 + x^2","baseline":"Abs[x^2 + 1]","note":"x^2+1>0 for real x (deep pass)","hyp":"pass"},
{"id":"abs15","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Sign[-x^2 - 1], Element[x, Reals]]","expected":"-1","baseline":"Sign[-x^2 - 1]","note":"-x^2-1<0 for real x (deep pass)","hyp":"pass"},
{"id":"abs16","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Sign[x^2 - x y + y^2 + 1], Element[x | y, Reals]]","expected":"1","baseline":"Sign[x^2 - x y + y^2 + 1]","note":"pos-definite quadratic form +1 (deep pass)","hyp":"pass"},
{"id":"abs17","cat":"2 Abs/Sign/Conj/Arg","expr":"Refine[Abs[Pi x], x > 0]","expected":"Pi x","baseline":"Abs[Pi x]","note":"Pi>0 and x>0 => |Pi x|=Pi x","hyp":"pass"},

# =====================================================================
# 3. Log / Exp branch behaviour
# =====================================================================
{"id":"log01","cat":"3 Log/Exp","expr":"Refine[Log[x], x < 0]","expected":"I Pi + Log[-x]","baseline":"Log[x]","note":"principal log of negative real","hyp":"pass"},
{"id":"log02","cat":"3 Log/Exp","expr":"Refine[Log[x^p], x > 0 && Element[p, Reals]]","expected":"p Log[x]","baseline":"Log[x^p]","note":"log(x^p)=p log x for x>0, real p","hyp":"pass"},
{"id":"log03","cat":"3 Log/Exp","expr":"Refine[Log[x^2], Element[x, Reals]]","expected":"2 Log[Abs[x]]","baseline":"Log[x^2]","note":"log(x^2)=2 log|x| for real x","hyp":"gap"},
{"id":"log04","cat":"3 Log/Exp","expr":"Refine[Exp[Log[x]], x > 0]","expected":"x","baseline":"Exp[Log[x]]","note":"e^(log x)=x for x>0","hyp":"pass"},
{"id":"log05","cat":"3 Log/Exp","expr":"Refine[Log[E^x], Element[x, Reals]]","expected":"x","baseline":"Log[E^x]","note":"log(e^x)=x for real x","hyp":"gap"},
{"id":"log06","cat":"3 Log/Exp","expr":"Refine[Log[a b], a > 0 && b > 0]","expected":"Log[a] + Log[b]","baseline":"Log[a b]","note":"log(ab)=log a+log b for a,b>0","hyp":"gap"},
{"id":"log07","cat":"3 Log/Exp","expr":"Refine[Log[1/x], x > 0]","expected":"-Log[x]","baseline":"Log[1/x]","note":"log(1/x)=-log x for x>0","hyp":"gap"},
{"id":"log08","cat":"3 Log/Exp","expr":"Refine[Log[x^p], x > 0]","expected":"p Log[x]","baseline":"Log[x^p]","note":"x>0 real base; holds for real p","hyp":"gap"},
{"id":"log09","cat":"3 Log/Exp","expr":"Refine[Log[Sqrt[x]], x > 0]","expected":"Log[x]/2","baseline":"Log[Sqrt[x]]","note":"log sqrt x = (1/2) log x for x>0","hyp":"gap"},

# =====================================================================
# 4. Integer / even / odd trig & parity
# =====================================================================
{"id":"trg01","cat":"4 Integer trig/parity","expr":"Refine[Sin[k Pi], Element[k, Integers]]","expected":"0","baseline":"Sin[k Pi]","note":"sin(k pi)=0","hyp":"pass"},
{"id":"trg02","cat":"4 Integer trig/parity","expr":"Refine[Cos[k Pi], Element[k, Integers]]","expected":"(-1)^k","baseline":"Cos[k Pi]","note":"cos(k pi)=(-1)^k","hyp":"pass"},
{"id":"trg03","cat":"4 Integer trig/parity","expr":"Refine[Tan[k Pi], Element[k, Integers]]","expected":"0","baseline":"Tan[k Pi]","note":"tan(k pi)=0","hyp":"pass"},
{"id":"trg04","cat":"4 Integer trig/parity","expr":"Refine[Cos[x + k Pi], Element[k, Integers]]","expected":"(-1)^k Cos[x]","baseline":"Cos[x + k Pi]","note":"cos shift by k pi","hyp":"pass"},
{"id":"trg05","cat":"4 Integer trig/parity","expr":"Refine[Sin[x + k Pi], Element[k, Integers]]","expected":"(-1)^k Sin[x]","baseline":"Sin[x + k Pi]","note":"sin shift by k pi","hyp":"pass"},
{"id":"trg06","cat":"4 Integer trig/parity","expr":"Refine[Tan[x + k Pi], Element[k, Integers]]","expected":"Tan[x]","baseline":"Tan[x + k Pi]","note":"tan has period pi","hyp":"pass"},
{"id":"trg07","cat":"4 Integer trig/parity","expr":"Refine[(-1)^(2 k), Element[k, Integers]]","expected":"1","baseline":"(-1)^(2 k)","note":"(-1)^(2k)=1","hyp":"pass"},
{"id":"trg08","cat":"4 Integer trig/parity","expr":"Refine[Exp[2 Pi I k], Element[k, Integers]]","expected":"1","baseline":"Exp[2 Pi I k]","note":"e^(2 pi i k)=1 for integer k","hyp":"gap"},
{"id":"trg09","cat":"4 Integer trig/parity","expr":"Refine[Sin[(2 k + 1) Pi/2], Element[k, Integers]]","expected":"(-1)^k","baseline":"Sin[(2 k + 1) Pi/2]","note":"sin((2k+1)pi/2)=(-1)^k","hyp":"gap"},
{"id":"trg10","cat":"4 Integer trig/parity","expr":"Refine[Cos[2 k Pi], Element[k, Integers]]","expected":"1","baseline":"Cos[2 k Pi]","note":"cos(2k pi)=1","hyp":"gap"},
{"id":"trg11","cat":"4 Integer trig/parity","expr":"Refine[ArcTan[Tan[x]], -Pi/2 < Re[x] < Pi/2]","expected":"x","baseline":"ArcTan[Tan[x]]","note":"arctan(tan x)=x on principal strip","hyp":"pass"},
{"id":"trg12","cat":"4 Integer trig/parity","expr":"Refine[(-1)^m, Element[m, Integers] && Mod[m, 2] == 0]","expected":"1","baseline":"(-1)^m","note":"even m => (-1)^m=1","hyp":"pass"},

# =====================================================================
# 5. Floor / Ceiling / Round / IntegerPart / FractionalPart / Mod
# =====================================================================
{"id":"flr01","cat":"5 Floor/Ceil/Mod","expr":"Refine[Floor[n], Element[n, Integers]]","expected":"n","baseline":"Floor[n]","note":"floor of integer","hyp":"pass"},
{"id":"flr02","cat":"5 Floor/Ceil/Mod","expr":"Refine[Floor[2 a + 1], Element[a, Integers]]","expected":"1 + 2 a","baseline":"Floor[2 a + 1]","note":"floor of integer expr","hyp":"pass"},
{"id":"flr03","cat":"5 Floor/Ceil/Mod","expr":"Refine[Ceiling[x], 2 < x <= 3]","expected":"3","baseline":"Ceiling[x]","note":"ceil constant =3 on (2,3]","hyp":"pass"},
{"id":"flr04","cat":"5 Floor/Ceil/Mod","expr":"Refine[Floor[x], 2 < x < 3]","expected":"2","baseline":"Floor[x]","note":"floor constant =2 on (2,3)","hyp":"pass"},
{"id":"flr05","cat":"5 Floor/Ceil/Mod","expr":"Refine[Ceiling[n], Element[n, Integers]]","expected":"n","baseline":"Ceiling[n]","note":"ceil of integer","hyp":"pass"},
{"id":"flr06","cat":"5 Floor/Ceil/Mod","expr":"Refine[Round[n], Element[n, Integers]]","expected":"n","baseline":"Round[n]","note":"round of integer","hyp":"pass"},
{"id":"flr07","cat":"5 Floor/Ceil/Mod","expr":"Refine[IntegerPart[n], Element[n, Integers]]","expected":"n","baseline":"IntegerPart[n]","note":"integerpart of integer","hyp":"pass"},
{"id":"flr08","cat":"5 Floor/Ceil/Mod","expr":"Refine[FractionalPart[n], Element[n, Integers]]","expected":"0","baseline":"FractionalPart[n]","note":"fracpart of integer=0","hyp":"pass"},
{"id":"flr09","cat":"5 Floor/Ceil/Mod","expr":"Refine[FractionalPart[a], a < 0 && Mod[a, 1] == 1/3]","expected":"-2/3","baseline":"FractionalPart[a]","note":"frac part from mod fact, a<0","hyp":"pass"},
{"id":"flr10","cat":"5 Floor/Ceil/Mod","expr":"Refine[Mod[a, 4], Element[(a + 3)/4, Integers]]","expected":"1","baseline":"Mod[a, 4]","note":"a=-3 mod 4 => a mod 4=1","hyp":"pass"},
{"id":"flr11","cat":"5 Floor/Ceil/Mod","expr":"Refine[Mod[n, 1], Element[n, Integers]]","expected":"0","baseline":"Mod[n, 1]","note":"integer mod 1 =0","hyp":"pass"},
{"id":"flr12","cat":"5 Floor/Ceil/Mod","expr":"Refine[Floor[x], 2 < x <= 3]","expected":"Floor[x]","baseline":"Floor[x]","note":"NOT constant on (2,3]: 2 on (2,3), 3 at 3 -> must stay (sound)","hyp":"sound"},
{"id":"flr13","cat":"5 Floor/Ceil/Mod","expr":"Refine[FractionalPart[2 n], Element[n, Integers]]","expected":"0","baseline":"FractionalPart[2 n]","note":"2n integer => fracpart=0","hyp":"pass"},

# =====================================================================
# 6. Re / Im / Conjugate / ComplexExpand
# =====================================================================
{"id":"cpx01","cat":"6 Re/Im/ComplexExpand","expr":"Refine[Re[a + b I], Element[a | b, Reals]]","expected":"a","baseline":"Re[a + b I]","note":"Re(a+bi)=a","hyp":"pass"},
{"id":"cpx02","cat":"6 Re/Im/ComplexExpand","expr":"Refine[Im[a + b I], Element[a | b, Reals]]","expected":"b","baseline":"Im[a + b I]","note":"Im(a+bi)=b","hyp":"pass"},
{"id":"cpx03","cat":"6 Re/Im/ComplexExpand","expr":"Refine[Re[x], Element[x, Reals]]","expected":"x","baseline":"Re[x]","note":"Re(real)=itself","hyp":"pass"},
{"id":"cpx04","cat":"6 Re/Im/ComplexExpand","expr":"Refine[Im[x], Element[x, Reals]]","expected":"0","baseline":"Im[x]","note":"Im(real)=0","hyp":"pass"},
{"id":"cpx05","cat":"6 Re/Im/ComplexExpand","expr":"Refine[Re[(a + b I)^2], Element[a | b, Reals]]","expected":"a^2 - b^2","baseline":"Re[(a + b I)^2]","note":"Re((a+bi)^2)=a^2-b^2","hyp":"pass"},
{"id":"cpx06","cat":"6 Re/Im/ComplexExpand","expr":"Refine[Im[(a + b I)^2], Element[a | b, Reals]]","expected":"2 a b","baseline":"Im[(a + b I)^2]","note":"Im((a+bi)^2)=2ab","hyp":"pass"},
{"id":"cpx07","cat":"6 Re/Im/ComplexExpand","expr":"Refine[Abs[a + b I], Element[a | b, Reals]]","expected":"Sqrt[a^2 + b^2]","baseline":"Abs[a + b I]","note":"|a+bi|=sqrt(a^2+b^2)","hyp":"gap"},
{"id":"cpx08","cat":"6 Re/Im/ComplexExpand","expr":"Refine[Conjugate[x^2 + I x], Element[x, Reals]]","expected":"x^2 - I x","baseline":"Conjugate[x^2 + I x]","note":"conj poly with real x","hyp":"pass"},
{"id":"cpx09","cat":"6 Re/Im/ComplexExpand","expr":"Refine[Re[a + b I], a > 0 && Element[b, Reals]]","expected":"a","baseline":"Re[a + b I]","note":"a>0 implies a real; Re=a","hyp":"pass"},
{"id":"cpx10","cat":"6 Re/Im/ComplexExpand","expr":"Refine[Im[a I], Element[a, Reals]]","expected":"a","baseline":"Im[a I]","note":"Im(a i)=a for real a","hyp":"pass"},

# =====================================================================
# 7. Element / domain decisions
# =====================================================================
{"id":"elt01","cat":"7 Element/domain","expr":"Refine[Element[k, Reals], Element[k, Integers]]","expected":"True","baseline":None,"note":"Integers subset Reals","hyp":"pass"},
{"id":"elt02","cat":"7 Element/domain","expr":"Refine[Element[k, Rationals], Element[k, Integers]]","expected":"True","baseline":None,"note":"Integers subset Rationals","hyp":"pass"},
{"id":"elt03","cat":"7 Element/domain","expr":"Refine[Element[k, Algebraics], Element[k, Integers]]","expected":"True","baseline":None,"note":"Integers subset Algebraics","hyp":"pass"},
{"id":"elt04","cat":"7 Element/domain","expr":"Refine[Element[k, Complexes], Element[k, Reals]]","expected":"True","baseline":None,"note":"Reals subset Complexes","hyp":"pass"},
{"id":"elt05","cat":"7 Element/domain","expr":"Refine[Element[x, Integers], Element[x, Reals]]","expected":"Element[x, Integers]","baseline":"Element[x, Integers]","note":"real !=> integer; undecidable stays (sound)","hyp":"sound"},
{"id":"elt06","cat":"7 Element/domain","expr":"Refine[Element[(2 x + x^p)/(x Gamma[x + 2]), Reals], x > 0 && p > 0]","expected":"True","baseline":None,"note":"real closed under the ops for x,p>0","hyp":"pass"},
{"id":"elt07","cat":"7 Element/domain","expr":"Refine[Element[2 k^3 Floor[x]^k, Integers], Element[k, Integers] && k > 0 && Element[x, Reals]]","expected":"True","baseline":None,"note":"integer combination","hyp":"pass"},
{"id":"elt08","cat":"7 Element/domain","expr":"Refine[Element[x, Reals], x^2 < 1]","expected":"True","baseline":None,"note":"algebraic-in-inequality => real","hyp":"pass"},
{"id":"elt09","cat":"7 Element/domain","expr":"Refine[Element[Sqrt[2], Algebraics], True]","expected":"True","baseline":"Element[Sqrt[2], Algebraics]","note":"sqrt2 is algebraic","hyp":"gap"},
{"id":"elt10","cat":"7 Element/domain","expr":"Refine[Element[x, Positive], x > 0]","expected":"True","baseline":"Element[x, Positive]","note":"Positive as queried domain","hyp":"gap"},
{"id":"elt11","cat":"7 Element/domain","expr":"Refine[Element[x, Negative], x < 0]","expected":"True","baseline":"Element[x, Negative]","note":"Negative as queried domain","hyp":"gap"},
{"id":"elt12","cat":"7 Element/domain","expr":"Refine[Element[x, NonNegative], x >= 0]","expected":"True","baseline":"Element[x, NonNegative]","note":"NonNegative as queried domain","hyp":"gap"},
{"id":"elt13","cat":"7 Element/domain","expr":"Refine[Element[n + m, Integers], Element[n, Integers] && Element[m, Integers]]","expected":"True","baseline":None,"note":"sum of integers","hyp":"pass"},
{"id":"elt14","cat":"7 Element/domain","expr":"Refine[Element[n m, Integers], Element[n | m, Integers]]","expected":"True","baseline":None,"note":"product of integers","hyp":"pass"},
{"id":"elt15","cat":"7 Element/domain","expr":"Refine[Element[7, Primes], True]","expected":"True","baseline":None,"note":"7 is prime (numeric)","hyp":"pass"},
{"id":"elt16","cat":"7 Element/domain","expr":"Refine[Element[k^2, NonNegative], Element[k, Reals]]","expected":"True","baseline":"Element[k^2, NonNegative]","note":"square of real is nonneg","hyp":"gap"},

# =====================================================================
# 8. Equal / Unequal predicate folding
# =====================================================================
{"id":"eqp01","cat":"8 Equal/Unequal","expr":"Refine[a^2 - b^2 + 1 == 0, a + b == 0]","expected":"False","baseline":None,"note":"a=-b => a^2-b^2=0 =>1=0 false","hyp":"pass"},
{"id":"eqp02","cat":"8 Equal/Unequal","expr":"Refine[x^2 + 1 == 0, Element[x, Reals]]","expected":"False","baseline":None,"note":"no real root","hyp":"pass"},
{"id":"eqp03","cat":"8 Equal/Unequal","expr":"Refine[x^2 == 1, x == 1]","expected":"True","baseline":None,"note":"x=1 => x^2=1","hyp":"pass"},
{"id":"eqp04","cat":"8 Equal/Unequal","expr":"Refine[a == b, a - b == 0]","expected":"True","baseline":None,"note":"a-b=0 => a=b","hyp":"pass"},
{"id":"eqp05","cat":"8 Equal/Unequal","expr":"Refine[a == b, a > b]","expected":"False","baseline":None,"note":"a>b => a!=b","hyp":"pass"},
{"id":"eqp06","cat":"8 Equal/Unequal","expr":"Refine[Sin[k Pi] == 0, Element[k, Integers]]","expected":"True","baseline":None,"note":"sin(k pi)=0 for integer k","hyp":"pass"},
{"id":"eqp07","cat":"8 Equal/Unequal","expr":"Refine[a b == 0, a == 0]","expected":"True","baseline":None,"note":"a=0 => ab=0","hyp":"pass"},
{"id":"eqp08","cat":"8 Equal/Unequal","expr":"Refine[x^2 != 1, x == 2]","expected":"True","baseline":None,"note":"x=2 => x^2=4 !=1","hyp":"pass"},
{"id":"eqp09","cat":"8 Equal/Unequal","expr":"Refine[x^2 - y^2 == (x - y) (x + y), True]","expected":"True","baseline":None,"note":"algebraic identity via zero test, no assumptions","hyp":"pass"},
{"id":"eqp10","cat":"8 Equal/Unequal","expr":"Refine[x == y, x^2 == y^2]","expected":"Refine[x == y, x^2 == y^2]","baseline":None,"note":"x^2=y^2 => x=+-y, cannot decide x=y; undecided (sound)","hyp":"sound"},

# =====================================================================
# 9. Inequalities & logic (CAD entailment)
# =====================================================================
{"id":"inq01","cat":"9 Inequalities/logic","expr":"Refine[a^2 - a b + b^2 >= 0, Element[a | b, Reals]]","expected":"True","baseline":None,"note":"pos-semidefinite form","hyp":"pass"},
{"id":"inq02","cat":"9 Inequalities/logic","expr":"Refine[(x - 1)^2 + (y - 2)^2 < 3/2, x^2 + y^2 <= 1]","expected":"False","baseline":None,"note":"disk containment contradiction","hyp":"pass"},
{"id":"inq03","cat":"9 Inequalities/logic","expr":"Refine[-1 < x < 1, x^2 < 1]","expected":"True","baseline":None,"note":"x^2<1 <=> -1<x<1","hyp":"pass"},
{"id":"inq04","cat":"9 Inequalities/logic","expr":"Refine[x^2 >= 0, Element[x, Reals]]","expected":"True","baseline":None,"note":"square nonneg","hyp":"pass"},
{"id":"inq05","cat":"9 Inequalities/logic","expr":"Refine[x^3 > 0, x > 0]","expected":"True","baseline":None,"note":"x>0 => x^3>0","hyp":"pass"},
{"id":"inq06","cat":"9 Inequalities/logic","expr":"Refine[x + y > 0, x > 0 && y > 0]","expected":"True","baseline":None,"note":"sum of positives","hyp":"pass"},
{"id":"inq07","cat":"9 Inequalities/logic","expr":"Refine[a > c, a > b && b > c]","expected":"True","baseline":None,"note":"transitivity","hyp":"pass"},
{"id":"inq08","cat":"9 Inequalities/logic","expr":"Refine[a >= c, a >= b && b >= c]","expected":"True","baseline":None,"note":"transitivity (non-strict)","hyp":"pass"},
{"id":"inq09","cat":"9 Inequalities/logic","expr":"Refine[x^2 + y^2 + z^2 >= 0, Element[x | y | z, Reals]]","expected":"True","baseline":None,"note":"sum of squares nonneg (3 vars)","hyp":"pass"},
{"id":"inq10","cat":"9 Inequalities/logic","expr":"Refine[Implies[x > 0, x^2 > 0], Element[x, Reals]]","expected":"True","baseline":None,"note":"valid implication","hyp":"pass"},
{"id":"inq11","cat":"9 Inequalities/logic","expr":"Refine[x > 0 || x <= 0, Element[x, Reals]]","expected":"True","baseline":None,"note":"tautology over reals","hyp":"pass"},
{"id":"inq12","cat":"9 Inequalities/logic","expr":"Refine[Not[x^2 < 0], Element[x, Reals]]","expected":"True","baseline":None,"note":"no real x with x^2<0","hyp":"pass"},
{"id":"inq13","cat":"9 Inequalities/logic","expr":"Refine[x y > 0, x > 0 && y > 0]","expected":"True","baseline":None,"note":"product of positives","hyp":"pass"},
{"id":"inq14","cat":"9 Inequalities/logic","expr":"Refine[x < 1, x < 0]","expected":"True","baseline":None,"note":"x<0 => x<1","hyp":"pass"},
{"id":"inq15","cat":"9 Inequalities/logic","expr":"Refine[x > 10, x > 5]","expected":"Refine[x > 10, x > 5]","baseline":None,"note":"x>5 neither entails nor refutes x>10; undecided (sound)","hyp":"sound"},
{"id":"inq16","cat":"9 Inequalities/logic","expr":"Refine[a + b + c + d + e + f + g > 0, a > 0 && b > 0 && c > 0 && d > 0 && e > 0 && f > 0 && g > 0]","expected":"True","baseline":"a + b + c + d + e + f + g > 0","note":"true, but 7 vars -> CAD bail (documents >6 limit)","hyp":"limit"},

# =====================================================================
# 10. Deep positivity (compound Sign/Abs/Sqrt)
# =====================================================================
{"id":"dpp01","cat":"10 Deep positivity","expr":"Refine[Sqrt[(x^2 + 1)^2], Element[x, Reals]]","expected":"1 + x^2","baseline":"Sqrt[(x^2 + 1)^2]","note":"sqrt of square of positive","hyp":"pass"},
{"id":"dpp02","cat":"10 Deep positivity","expr":"Refine[Sign[x^2 + 1], Element[x, Reals]]","expected":"1","baseline":"Sign[x^2 + 1]","note":"REAL x: x^2+1>0 (deep pass). Without a real assumption x could be complex (x=I gives 0), so the no-assumption form correctly stays; this case now supplies the real fact","hyp":"pass"},
{"id":"dpp03","cat":"10 Deep positivity","expr":"Refine[Abs[x^2 - 2 x + 1], Element[x, Reals]]","expected":"1 - 2 x + x^2","baseline":"Abs[x^2 - 2 x + 1]","note":"(x-1)^2>=0 (nonneg, not strict) -> deep pass only does strict","hyp":"gap"},
{"id":"dpp04","cat":"10 Deep positivity","expr":"Refine[Sqrt[(x - 1)^2], Element[x, Reals]]","expected":"Abs[-1 + x]","baseline":"Sqrt[(x - 1)^2]","note":"sqrt((x-1)^2)=|x-1| for real x (compound base)","hyp":"gap"},
{"id":"dpp05","cat":"10 Deep positivity","expr":"Refine[Sign[x^3], x > 0]","expected":"1","baseline":"Sign[x^3]","note":"x>0 => x^3>0","hyp":"pass"},
{"id":"dpp06","cat":"10 Deep positivity","expr":"Refine[Abs[x^4 + x^2 + 1], Element[x, Reals]]","expected":"1 + x^2 + x^4","baseline":"Abs[x^4 + x^2 + 1]","note":"strictly positive quartic","hyp":"pass"},
{"id":"dpp07","cat":"10 Deep positivity","expr":"Refine[Sign[(x - 2) (x - 3)], x > 3]","expected":"1","baseline":"Sign[(x - 2) (x - 3)]","note":"x>3 => both factors >0","hyp":"pass"},
{"id":"dpp08","cat":"10 Deep positivity","expr":"Refine[Abs[2 x + 1], x > 0]","expected":"1 + 2 x","baseline":"Abs[2 x + 1]","note":"2x+1>0 for x>0 (compound)","hyp":"pass"},

# =====================================================================
# 11. Assumption plumbing ($Assumptions, Assuming, options, lists)
# =====================================================================
{"id":"plm01","cat":"11 Plumbing","expr":"Assuming[x > 0, Refine[Sqrt[x^2]]]","expected":"x","baseline":"Sqrt[x^2]","note":"Assuming sets $Assumptions","hyp":"pass"},
{"id":"plm02","cat":"11 Plumbing","expr":"Assuming[x > 0, Refine[Sqrt[x^2 y^2], y < 0]]","expected":"-x y","baseline":None,"note":"positional combined with $Assumptions","hyp":"pass"},
{"id":"plm03","cat":"11 Plumbing","expr":"Assuming[x > 0, Refine[Sqrt[x^2 y^2], Assumptions -> y < 0]]","expected":"-Sqrt[x^2] y","baseline":None,"note":"option overrides $Assumptions: only y<0 known; x not even real so Sqrt[x^2] must stay (corpus expectation corrected from -Abs[x] y, which wrongly assumed x real)","hyp":"pass"},
{"id":"plm04","cat":"11 Plumbing","expr":"Refine[Cos[k Pi]^m, Element[k, Integers], Assumptions -> Mod[m, 2] == 0]","expected":"1","baseline":None,"note":"positional + option combine; (-1)^(k m), m even","hyp":"pass"},
{"id":"plm05","cat":"11 Plumbing","expr":"Refine[Sqrt[x^2]]","expected":"Sqrt[x^2]","baseline":"Sqrt[x^2]","note":"no assumptions -> identity","hyp":"pass"},
{"id":"plm06","cat":"11 Plumbing","expr":"Refine[Sqrt[x^2], True]","expected":"Sqrt[x^2]","baseline":"Sqrt[x^2]","note":"True -> identity","hyp":"pass"},
{"id":"plm07","cat":"11 Plumbing","expr":"Refine[Sqrt[x^2], {x > 0}]","expected":"x","baseline":"Sqrt[x^2]","note":"list assumption","hyp":"pass"},
{"id":"plm08","cat":"11 Plumbing","expr":"Refine[Sqrt[x^2 y^2], {x > 0, y < 0}]","expected":"-x y","baseline":"Sqrt[x^2 y^2]","note":"list of two assumptions","hyp":"pass"},
{"id":"plm09","cat":"11 Plumbing","expr":"Refine[Sqrt[x^2], And[x > 0]]","expected":"x","baseline":"Sqrt[x^2]","note":"And wrapper","hyp":"pass"},
{"id":"plm10","cat":"11 Plumbing","expr":"Block[{$Assumptions = x > 0}, Refine[Sqrt[x^2]]]","expected":"x","baseline":"Sqrt[x^2]","note":"global $Assumptions read","hyp":"pass"},
{"id":"plm11","cat":"11 Plumbing","expr":"Refine[Sqrt[x^2 y^2 z^2], x > 0 && y < 0 && z > 0]","expected":"-x y z","baseline":"Sqrt[x^2 y^2 z^2]","note":"three sign facts","hyp":"pass"},
{"id":"plm12","cat":"11 Plumbing","expr":"Assuming[x > 0, Assuming[y < 0, Refine[Sqrt[x^2 y^2]]]]","expected":"-x y","baseline":None,"note":"nested Assuming compose","hyp":"pass"},
{"id":"plm13","cat":"11 Plumbing","expr":"Refine[Sqrt[x^2], Element[x, Reals] && x != 0]","expected":"Abs[x]","baseline":"Sqrt[x^2]","note":"real + Unequal fact","hyp":"pass"},
{"id":"plm14","cat":"11 Plumbing","expr":"Refine[Sqrt[x^2], x > 0, TimeConstraint -> 5]","expected":"x","baseline":"Sqrt[x^2]","note":"TimeConstraint option not positional","hyp":"pass"},

# =====================================================================
# 12. Adversarial / soundness / robustness
# =====================================================================
{"id":"adv01","cat":"12 Adversarial/soundness","expr":"Refine[Sqrt[x^2], False]","expected":"Sqrt[x^2]","baseline":"Sqrt[x^2]","note":"lone False dropped by ctx_walk -> identity (ctx->inconsistent unused)","hyp":"sound"},
{"id":"adv02","cat":"12 Adversarial/soundness","expr":"Refine[Abs[x], x > 0 && x < 0]","expected":"x","baseline":"Abs[x]","note":"contradictory premise; currently uses x>0 -> x (inconsistency undetected)","hyp":"sound"},
{"id":"adv03","cat":"12 Adversarial/soundness","expr":"Refine[Sign[(Sqrt[2] + Sqrt[3])^2 - 5 - 2 Sqrt[6]], True]","expected":"0","baseline":"Sign[(Sqrt[2] + Sqrt[3])^2 - 5 - 2 Sqrt[6]]","note":"argument is EXACTLY 0; ideal 0, must NOT be +-1 (soundness)","hyp":"sound"},
{"id":"adv04","cat":"12 Adversarial/soundness","expr":"Refine[Abs[x], x > 9007199254740993]","expected":"x","baseline":"Abs[x]","note":"bound > 2^53 (double precision)","hyp":"pass"},
{"id":"adv05","cat":"12 Adversarial/soundness","expr":"Refine[Abs[x], x < -9007199254740993]","expected":"-x","baseline":"Abs[x]","note":"large negative bound","hyp":"pass"},
{"id":"adv06","cat":"12 Adversarial/soundness","expr":"Refine[Abs[a1], a1 > 0 && a2 > 0 && a3 > 0 && a4 > 0 && a5 > 0 && a6 > 0 && a7 > 0 && a8 > 0 && a9 > 0 && a10 > 0 && a11 > 0 && a12 > 0 && a13 > 0 && a14 > 0 && a15 > 0 && a16 > 0 && a17 > 0 && a18 > 0 && a19 > 0 && a20 > 0]","expected":"a1","baseline":"Abs[a1]","note":"20 positive symbols: MAX_SYM=16 / 8192-byte buffer overflow drops rules?","hyp":"gap"},
{"id":"adv07","cat":"12 Adversarial/soundness","expr":"Refine[x, 5]","expected":"x","baseline":"x","note":"garbage numeric assumption -> identity, no crash","hyp":"pass"},
{"id":"adv08","cat":"12 Adversarial/soundness","expr":"Refine[x, \"foo\"]","expected":"x","baseline":"x","note":"string assumption -> identity, no crash","hyp":"pass"},
{"id":"adv09","cat":"12 Adversarial/soundness","expr":"Refine[]","expected":"Refine[]","baseline":None,"note":"0 args -> message + unevaluated","hyp":"pass"},
{"id":"adv10","cat":"12 Adversarial/soundness","expr":"Refine[a, b, c]","expected":"Refine[a, b, c]","baseline":None,"note":">2 args -> message + unevaluated","hyp":"pass"},
{"id":"adv11","cat":"12 Adversarial/soundness","expr":"Refine[x^2 + y^2 + z^2 + u^2 + v^2 + w^2 >= 0, Element[x | y | z | u | v | w, Reals]]","expected":"True","baseline":None,"note":"6-var sum of squares (at the CAD var ceiling)","hyp":"pass"},
{"id":"adv12","cat":"12 Adversarial/soundness","expr":"Refine[x^2 + y^2 + z^2 + u^2 + v^2 + w^2 >= 0, Element[x | y | z | u | v | w, Reals] && TimeConstraint -> 1]","expected":"True","baseline":None,"note":"placeholder; real TC case below","hyp":"pass"},
{"id":"adv13","cat":"12 Adversarial/soundness","expr":"Refine[Sign[0], True]","expected":"0","baseline":None,"note":"sign(0)=0","hyp":"pass"},
{"id":"adv14","cat":"12 Adversarial/soundness","expr":"Refine[Abs[-5], True]","expected":"5","baseline":None,"note":"numeric passthrough","hyp":"pass"},
{"id":"adv15","cat":"12 Adversarial/soundness","expr":"Refine[Nest[Sqrt[#^2] &, x, 3], x > 0]","expected":"x","baseline":None,"note":"nested sqrt-of-square, x>0","hyp":"gap"},
{"id":"adv16","cat":"12 Adversarial/soundness","expr":"Refine[Sqrt[x^2], x > 0, TimeConstraint -> Infinity]","expected":"x","baseline":"Sqrt[x^2]","note":"TimeConstraint Infinity accepted","hyp":"pass"},
]
