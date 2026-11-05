(* RoachCorpus.m
 *
 * The test suites of Kelly Roach's three papers on hypergeometric and Meijer G
 * methods, as one scored corpus (runner: tests/test_roach_corpus.c):
 *
 *   R08-01 .. R08-30   "Solving Integrals with the Quantum Computer Algebra
 *                      System" (2008): the 30 definite integrals (pages 2-5:
 *                      two in the introduction, the 22 of the two tables, the
 *                      6 elliptic integrals);
 *   R08-F01 .. R08-F11 the same paper's 11 pFq special values (page 6);
 *   I96-01 .. I96-12   "Hypergeometric Function Representations" (ISSAC '96),
 *                      the sec. 21 gallery of pFq representations;
 *   I97-01 .. I97-12   "Meijer G Function Representations" (ISSAC '97), the
 *                      sec. 15 gallery of definite integrals (PBM references).
 *
 * Each record is
 *
 *   {id, Hold[expr], {param -> value, ...}, reference, {uneval heads}}
 *
 * expr is the question exactly as the paper poses it (an Integrate with the
 * assumptions its formula needs, or FunctionExpand of the pFq); the rules pick
 * a generic sample point inside those assumptions; reference is the value of
 * expr at that point to 30 digits; a result that still contains one of the
 * uneval heads has not been answered.
 *
 * THE REFERENCES ARE INDEPENDENT OF MATHILDA. Every one is computed by
 * tools/roach_corpus_refs.py with mpmath at 40 digits from the definition --
 * adaptive quadrature of the integrand, on a contour rotated into a half plane
 * of decay where the real-axis integral is only conditionally convergent, and
 * mpmath.hyper summing the pFq series. The paper's closed form is evaluated
 * there as a cross-check only; it agrees to 15+ digits on 59 of the 65 records.
 * The six where it does not, every one settled in favour of the quadrature:
 *   - R08-01  the printed form with the paper's Pi(n, k), n = 25/21 > 1, is
 *             complex in mpmath's convention; its real part (the principal
 *             value) agrees.
 *   - R08-24  Ei(-x^2) e^(a^2 x^2): the printed pi^(3/2)/(2a) - Sqrt[pi]
 *             ArcSin[a]/a is positive, the integrand negative; the value is
 *             -Sqrt[pi] ArcSin[a]/a (-> -Sqrt[pi] = int Ei(-x^2) as a -> 0).
 *   - R08-30  the printed F/E/Pi combination misses by 8% (by 1% with the
 *             other sign of the Pi characteristic); two quadratures agree to
 *             20 digits.
 *   - I96-10  3F2(1,2,3; 1/2,4; z): the first printed term has (1 - z) where
 *             the series needs (z - 1).
 *   - I97-06  x^(3/2) J0(b x) J0(c x) diverges as an ordinary integral; the
 *             reference is its Abel/Mellin value, the analytic continuation of
 *             the Weber-Schafheitlin 2F1 in the power of x (confirmed by the
 *             e^(-eps x) regularised integrals, -0.092, -0.120, -0.133 at
 *             eps = 0.2, 0.1, 0.05, extrapolating to -0.146); the printed
 *             elliptic form gives -0.127.
 *   - R08-20  (no disagreement) ArcTan[x]/(1 - x^2) has a pole at x = 1: the
 *             record asks for the principal value.
 * The cross-check also caught a reference of our own: one tanh-sinh pass over
 * (0, oo) for Sin[b x^2] K_nu(c x) (I97-11, and I97-10) is wrong at 2e-7; the interval is
 * now cut every 1/2, and the value agrees with the paper to 16 digits.
 * Keep the sample points in step with tools/roach_corpus_refs.py.
 *)

{
  (* ---- Roach 2008, the introduction ---- *)
  {"R08-01", Hold[Integrate[((z^2 + 1) (z^2 - 4))^(-1/2)/(z - 5), {z, -Infinity, -3}]], {}, -0.0301578041634659035640552445734, {"Integrate"}},
  {"R08-02", Hold[Integrate[E^(-al x) BesselJ[0, be x]^2, {x, 0, Infinity}, Assumptions -> al > 0 && be > 0]], {al -> 3/2, be -> 1}, 0.508099680048529116034736961805, {"Integrate"}},

  (* ---- Roach 2008, table 1 (page 3) ---- *)
  {"R08-03", Hold[Integrate[(a^2 - x^2)^(n - 1/2), {x, 0, a}, Assumptions -> a > 0 && n > 0]], {a -> 2, n -> 1/3}, 1.77828809261779239153078252454, {"Integrate"}},
  {"R08-04", Hold[Integrate[E^(-b x - x^2/(4 a)), {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]], {a -> 2, b -> 1/2}, 1.31135908483759694308774246146, {"Integrate"}},
  {"R08-05", Hold[Integrate[x/(E^(3 x) - 1)^(1/3), {x, 0, Infinity}]], {}, 1.02976166065275577729043151696, {"Integrate"}},
  {"R08-06", Hold[Integrate[x E^(-x n)/(E^x + 1), {x, 0, Infinity}, Assumptions -> n > 0]], {n -> 1/2}, 0.33613762329112393978158594027, {"Integrate"}},
  {"R08-07", Hold[Integrate[ArcSin[x] x^(2 n - 1), {x, 0, 1}, Assumptions -> n > 0]], {n -> 3/4}, 0.464518094687237788273138996673, {"Integrate"}},
  {"R08-08", Hold[Integrate[x^(a - 1)/Sqrt[x^b + 1], {x, 0, Infinity}, Assumptions -> a > 0 && b > 2 a]], {a -> 3/4, b -> 4}, 1.98695810057561096269351789814, {"Integrate"}},
  {"R08-09", Hold[Integrate[Cos[x]^(-2 a - 2) Cos[2 x]^a, {x, 0, Pi/4}, Assumptions -> a > 0]], {a -> 1/3}, 0.84130926319527255670501144743, {"Integrate"}},
  {"R08-10", Hold[Integrate[x E^(-2 x n)/(E^x + 1), {x, 0, Infinity}, Assumptions -> n > 0]], {n -> 1/2}, 0.177532966575886781763792416677, {"Integrate"}},
  {"R08-11", Hold[Integrate[x E^(-3 x)/(E^(-x) + 1), {x, 0, Infinity}]], {}, 0.072467033424113218236207583323, {"Integrate"}},
  {"R08-12", Hold[Integrate[x^(b - 1) (a x + 1)^(-c), {x, 0, Infinity}, Assumptions -> a > 0 && b > 0 && c > b]], {a -> 2, b -> 1/2, c -> 3/2}, 1.41421356237309504880145495805, {"Integrate"}},
  {"R08-13", Hold[Integrate[(SinIntegral[a x] - Pi/2) E^(-b x), {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]], {a -> 2, b -> 1}, -0.463647609000806116214256231461, {"Integrate"}},

  (* ---- Roach 2008, table 2 (page 4) ---- *)
  {"R08-14", Hold[Integrate[x E^(-x)/Sqrt[E^x - 1], {x, 0, Infinity}]], {}, 0.606789763508705511269367206598, {"Integrate"}},
  {"R08-15", Hold[Integrate[x^(a - 1) (b + x)^(-n), {x, 0, Infinity}, Assumptions -> a > 0 && n > a && b > 0]], {a -> 1/2, n -> 2, b -> 3}, 0.302299894039036308432330789663, {"Integrate"}},
  {"R08-16", Hold[Integrate[BesselJ[1, a x] E^(-x^2 b), {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]], {a -> 1, b -> 1/2}, 0.393469340287366576396200465009, {"Integrate"}},
  {"R08-17", Hold[Integrate[E^(-a x) x^n/(b + x), {x, 0, Infinity}, Assumptions -> a > 0 && b > 0 && n > -1]], {a -> 1, b -> 2, n -> 1/2}, 0.278738824696113026082434476756, {"Integrate"}},
  {"R08-18", Hold[Integrate[ArcCot[x]/(x + 1), {x, 0, Infinity}]], {}, 1.46036211675311954767977573949, {"Integrate"}},
  {"R08-19", Hold[Integrate[Log[x + 1] Cos[a x]/x, {x, 0, Infinity}, Assumptions -> a > 0]], {a -> 1}, 0.25205402997274035110102394302, {"Integrate"}},
  {"R08-20", Hold[Integrate[ArcTan[x]/(1 - x^2), {x, 0, Infinity}, PrincipalValue -> True]], {}, -0.915965594177219015054603514932, {"Integrate"}},
  {"R08-21", Hold[Integrate[BesselJ[mu, a x]/(b + x), {x, 0, Infinity}, Assumptions -> a > 0 && b > 0 && mu > -1]], {mu -> 1/3, a -> 1, b -> 2}, 0.390476209707218492268201503129, {"Integrate"}},
  {"R08-22", Hold[Integrate[x^2 Csch[x] E^(-2 x n), {x, 0, Infinity}, Assumptions -> n > 0]], {n -> 1/4}, 1.32565830408440990785817672761, {"Integrate"}},
  {"R08-23", Hold[Integrate[x Erfc[a x] E^(-x^2 b), {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]], {a -> 1, b -> 1/2}, 0.183503419072273967267571975098, {"Integrate"}},
  {"R08-24", Hold[Integrate[ExpIntegralEi[-x^2] E^(x^2 a^2), {x, 0, Infinity}, Assumptions -> 0 < a < 1]], {a -> 1/2}, -1.85610933227723594842827266071, {"Integrate"}},

  (* ---- Roach 2008, the elliptic integrals (page 5) ---- *)
  {"R08-25", Hold[Integrate[1/((z + 14) Sqrt[(z - 1) ((z - 13)^2 + 81)]), {z, 2, Infinity}]], {}, 0.0355229827878928709690543991018, {"Integrate"}},
  {"R08-26", Hold[Integrate[1/((18 - z) Sqrt[(1 - z) ((16 - z)^2 + 64)]), {z, -Infinity, 0}]], {}, 0.0160745346391821514209538292326, {"Integrate"}},
  {"R08-27", Hold[Integrate[Sqrt[(z^2 - 1)/(9 - z^2)], {z, 1, 2}]], {}, 0.429353457153347287735897684271, {"Integrate"}},
  {"R08-28", Hold[Integrate[1/((z + 5) Sqrt[(z - 1) (5 - z) ((z - 1)^2 + 9)]), {z, 4, 5}]], {}, 0.0230406134564657240400508717783, {"Integrate"}},
  {"R08-29", Hold[Integrate[Sqrt[z/((z - 1) (z - 2))], {z, 2, 3}]], {}, 2.66963010670665880951270357469, {"Integrate"}},
  {"R08-30", Hold[Integrate[1/((6 - z) Sqrt[(5 - z) (z - 2) (z - 3) z^3]), {z, 3, 4}]], {}, 0.0855294314210272527170793977095, {"Integrate"}},

  (* ---- Roach 2008, the pFq special values (page 6); z -> 3/10 ---- *)
  {"R08-F01", Hold[FunctionExpand[HypergeometricPFQ[{-26/7, -26/7}, {1}, -1]]], {}, 4.58057136006074500144181536408, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"R08-F02", Hold[FunctionExpand[HypergeometricPFQ[{-15/4, 1, 1}, {-7/2, 2}, 1]]], {}, 5.53125475132864724067456977782, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"R08-F03", Hold[FunctionExpand[HypergeometricPFQ[{1, 1, 3/2}, {1/2, 9/2}, -1]]], {}, 0.581186265848152708105397604222, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"R08-F04", Hold[FunctionExpand[HypergeometricPFQ[{-1/2, -1/2, 5/2}, {3/2, 7/2}, -1]]], {}, 0.884944864116929866214372538296, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"R08-F05", Hold[FunctionExpand[HypergeometricPFQ[{-11/4, -5/2, -5/2}, {-7/2, 1/4}, 1]]], {}, 33.2439439919179475969649305147, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"R08-F06", Hold[FunctionExpand[HypergeometricPFQ[{1/2, 1/2}, {1}, 9/10]]], {}, 1.64126441434237073328699974737, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"R08-F07", Hold[FunctionExpand[HypergeometricPFQ[{1, 1}, {2/3}, 8/9]]], {}, 22.028357412103631581160968608, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"R08-F08", Hold[FunctionExpand[HypergeometricPFQ[{1, 1}, {1/2}, -7/8]]], {}, 0.229127500629833605910668080713, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"R08-F09", Hold[FunctionExpand[HypergeometricPFQ[{1, 3}, {2, 4, 5}, -z]]], {z -> 3/10}, 0.977797338564571728172753593452, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"R08-F10", Hold[FunctionExpand[HypergeometricPFQ[{-7/2}, {-5/2}, -z]]], {z -> 3/10}, 0.651183374660893249025833420772, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"R08-F11", Hold[FunctionExpand[HypergeometricPFQ[{2}, {9/4}, z]]], {z -> 3/10}, 1.30732809668195217742669362108, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},

  (* ---- ISSAC '96, sec. 21 gallery; z -> 3/10 ---- *)
  {"I96-01", Hold[FunctionExpand[HypergeometricPFQ[{1/2}, {9/2}, z]]], {z -> 3/10}, 1.03475131295689413470102504821, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"I96-02", Hold[FunctionExpand[HypergeometricPFQ[{-3/2}, {-1/2, 1/2}, z]]], {z -> 3/10}, 2.61756884689041260242559973663, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"I96-03", Hold[FunctionExpand[HypergeometricPFQ[{-3/2, -1/2}, {-5/2, 1}, z]]], {z -> 3/10}, 0.908939744244027533438705399424, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"I96-04", Hold[FunctionExpand[HypergeometricPFQ[{-3/2, -1/2}, {2}, z]]], {z -> 3/10}, 1.11106575300141736505769902197, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"I96-05", Hold[FunctionExpand[HypergeometricPFQ[{-1/2, 1}, {1/4, 1/2, 3/4}, z]]], {z -> 3/10}, -0.675308497035299942523285538048, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"I96-06", Hold[FunctionExpand[HypergeometricPFQ[{}, {1, 5/2, 4}, z]]], {z -> 3/10}, 1.03012873023605164896847427714, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"I96-07", Hold[FunctionExpand[HypergeometricPFQ[{3/2}, {5/2, 5}, z]]], {z -> 3/10}, 1.03665005509937087575890130669, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"I96-08", Hold[FunctionExpand[HypergeometricPFQ[{-1/2, 1, 2}, {3, 4}, z]]], {z -> 3/10}, 0.974400345321917486115325562032, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"I96-09", Hold[FunctionExpand[HypergeometricPFQ[{-1/2, 1/2, 1}, {3/2, 5/2}, z]]], {z -> 3/10}, 0.979444473960257876589927351008, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"I96-10", Hold[FunctionExpand[HypergeometricPFQ[{1, 2, 3}, {1/2, 4}, z]]], {z -> 3/10}, 2.60175080453983266532071738343, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"I96-11", Hold[FunctionExpand[HypergeometricPFQ[{1/3}, {-2/3, -1/2, 1/2, 1}, z]]], {z -> 3/10}, 1.84374477322918164993280239837, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},
  {"I96-12", Hold[FunctionExpand[HypergeometricPFQ[{1/2}, {-1/3, -1/2, -2/3}, z]]], {z -> 3/10}, -3.29226305627307202584447775507, {"HypergeometricPFQ", "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1", "MeijerG"}},

  (* ---- ISSAC '97, sec. 15 gallery (PBM numbers as printed) ---- *)
  (* 2.5.6(3) v1 *)
  {"I97-01", Hold[Integrate[Sin[b x]/(x^2 + c^2)^rho, {x, 0, Infinity}, Assumptions -> b > 0 && c > 0 && rho > 0]], {b -> 1, c -> 2, rho -> 3/4}, 0.376101470454130658755729023318, {"Integrate"}},
  (* 2.5.21(3b) v1 *)
  {"I97-02", Hold[Integrate[Cos[a x^2 + 2 b x], {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]], {a -> 1, b -> 1/2}, 0.270513580162214144258900856158, {"Integrate"}},
  (* 2.7.6(6) v1 *)
  {"I97-03", Hold[Integrate[Cos[b x] ArcTan[a/x^2], {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]], {a -> 1, b -> 1}, 1.00630124711007787401371982295, {"Integrate"}},
  (* 2.12.19(6) v2 *)
  {"I97-04", Hold[Integrate[Cos[b Sqrt[x]] BesselJ[0, c x]/Sqrt[x], {x, 0, Infinity}, Assumptions -> b > 0 && c > 0]], {b -> 1, c -> 1}, 1.83288316592168606882698947557, {"Integrate"}},
  (* 2.12.31(1l) v2 *)
  {"I97-05", Hold[Integrate[BesselJ[1, b x] BesselJ[1, c x]/x^2, {x, 0, Infinity}, Assumptions -> 0 < b < c]], {b -> 1, c -> 2}, 0.483843755630125818515495010646, {"Integrate"}},
  (* 2.12.31(1u) v2: Abel/Mellin value, see the header *)
  {"I97-06", Hold[Integrate[x^(3/2) BesselJ[0, b x] BesselJ[0, c x], {x, 0, Infinity}, Assumptions -> b > c > 0]], {b -> 2, c -> 1}, -0.144659295464170468725658691769, {"Integrate"}},
  (* 2.14.1(6a) v2 *)
  {"I97-07", Hold[Integrate[E^(I p x) HankelH1[0, c x], {x, 0, Infinity}, Assumptions -> 0 < p < c]], {p -> 1/2, c -> 1}, 0.769800358919501019345531707336, {"Integrate"}},
  (* 2.15.20(4c) v2 *)
  {"I97-08", Hold[Integrate[E^(-p x) BesselI[1, c x]^2, {x, 0, Infinity}, Assumptions -> c > 0 && p > 2 c]], {c -> 1, p -> 3}, 0.0280894710101649336251615385625, {"Integrate"}},
  (* 2.15.20(5f) v2 *)
  {"I97-09", Hold[Integrate[E^(-p x) BesselI[1, c x]^2/x, {x, 0, Infinity}, Assumptions -> c > 0 && p > 2 c]], {c -> 1, p -> 3}, 0.0339008480035941020747762752237, {"Integrate"}},
  (* 2.16.15(1a) v2 *)
  {"I97-10", Hold[Integrate[x^(nu + 1) Sin[c x^2/(2 a)] BesselK[nu, c x], {x, 0, Infinity}, Assumptions -> a > 0 && c > 0 && 0 < nu < 1]], {nu -> 1/3, a -> 1, c -> 1}, 0.373814491636919350022771522621, {"Integrate"}},
  (* 2.16.15(2a) v2 *)
  {"I97-11", Hold[Integrate[Sin[b x^2] BesselK[nu, c x], {x, 0, Infinity}, Assumptions -> b > 0 && c > 0 && 0 < nu < 1]], {nu -> 1/3, b -> 1, c -> 1}, 0.356536333834427244345020784515, {"Integrate"}},
  (* 2.7.16(3) v3 *)
  {"I97-12", Hold[Integrate[x BesselJ[-nu, b x] (BesselY[nu, c x] - StruveH[nu, c x]), {x, 0, Infinity}, Assumptions -> 0 < b < c && 0 < nu < 1/2]], {nu -> 1/4, b -> 1, c -> 2}, -0.178443761487844969414640840918, {"Integrate"}}
}
