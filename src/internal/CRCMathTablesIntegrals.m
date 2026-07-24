
(* Source: CRC Standard Mathematical Tables, 31st Edition *)
(*
   Loaded lazily on first call to Integrate's CRCTable stage (see
   src/integrate.c's try_crctable / crc_lazy_load).  The table is
   addressable by users as Integrate`CRCTable[f, x]; the internal rules
   are stored on the short name IntegrateTable for readability — the
   public wrapper at the bottom of this file forwards to it.

   Mathilda's BeginPackage/Begin parsing inside Get is incomplete (it
   mishandles SetDelayed inside an explicit context), so we keep the
   layout flat instead of using BeginPackage.
*)

(* Elementary forms *)
IntegrateTable[a_, x_] /; FreeQ[a, x] := a x;
IntegrateTable[a_ f_, x_] /; FreeQ[a, x] := a IntegrateTable[f, x];
IntegrateTable[f_ + g_, x_] := IntegrateTable[f, x] + IntegrateTable[g, x];
IntegrateTable[x_^n_., x_] /; FreeQ[n, x] && n =!= -1 := x^(n + 1)/(n + 1);
IntegrateTable[1/x_, x_] := Log[x];
IntegrateTable[Exp[a_. x_], x_] /; FreeQ[a, x] := Exp[a x]/a;
IntegrateTable[b_^(a_. x_), x_] /; FreeQ[{a, b}, x] && b > 0 := b^(a x)/(a Log[b]);
IntegrateTable[Log[a_. x_], x_] /; FreeQ[a, x] := x Log[a x] - x;
IntegrateTable[1/(x_^2 + a_), x_] /; FreeQ[a, x] && a > 0 := 1/Sqrt[a] ArcTan[x/Sqrt[a]]; 
IntegrateTable[1/(a_ - x_^2), x_] /; FreeQ[a, x] && a > 0 := 1/Sqrt[a] ArcTanh[x/Sqrt[a]];
IntegrateTable[1/(x_^2 + a_), x_] /; FreeQ[a, x] && a < 0 := -1/Sqrt[-a] ArcCoth[x/Sqrt[-a]];
IntegrateTable[Power[a_ - x_^2, -1/2], x_] /; FreeQ[a, x] && a > 0 := ArcSin[x/Sqrt[a]];
IntegrateTable[Power[x_^2 + a_, -1/2], x_] /; FreeQ[a, x] := Log[x + Sqrt[x^2 + a]];
IntegrateTable[Power[x_^2 - a_, -1/2], x_] /; FreeQ[a, x] := Log[x + Sqrt[x^2 - a]];
IntegrateTable[1/x_ Power[x_^2 - a_, -1/2], x_] /; FreeQ[a, x] := ArcTan[Sqrt[x^2 - a]/Sqrt[a]]/Sqrt[a];
IntegrateTable[1/x_ Power[b_. x_^2 + a_, -1/2], x_] /; FreeQ[{a, b}, x] := -ArcTanh[Sqrt[b x^2 + a]/Sqrt[a]]/Sqrt[a];

(* Forms containing a + b x *)
(* Formula 23 *)
IntegrateTable[(b_. x_ + a_.)^n_, x_] /; FreeQ[{a, b, n}, x] && n =!= -1 := ((a + b x)^(n + 1))/(b*(n + 1));

(* Formulae 24 *)
IntegrateTable[x_ (b_. x_ + a_.)^n_, x_] /; FreeQ[{a, b, n}, x] && n =!= -1 && n =!= -2 := (b x + a)^(n + 2)/(b^2 (n + 2)) - a (a + b x)^(n + 1)/(b^2 (n + 1));

(* Formula 25 *)
IntegrateTable[x_^2 (b_. x_ + a_.)^n_, x_] /; FreeQ[{a, b, n}, x] && n =!= -1 && n =!= -2 && n =!= -3 := 
  1/b^3 ((a + b x)^(n + 3)/(n + 3) - 2 a (a + b x)^(n + 2)/(n + 2) + a^2 (a + b x)^(n + 1)/(n + 1));

(* Formula 26: Three distinct reduction forms are provided. *)
(* Form 1: Reduces the power of the binomial *)
IntegrateTable[x_^m_ (b_. x_ + a_.)^n_, x_] /; FreeQ[{a, b, m, n}, x] && m + n + 1 =!= 0 && IntegerQ[m] && IntegerQ[n] && n > 0 :=
  (x^(m + 1) (a + b x)^n)/(m + n + 1) + (a n)/(m + n + 1) IntegrateTable[x^m (a + b x)^(n - 1), x];

(* Form 2: Increases the power of the binomial *)
IntegrateTable[x_^m_ (b_. x_ + a_.)^n_, x_] /; FreeQ[{a, b, m, n}, x] && n =!= -1 && IntegerQ[m] && IntegerQ[n] && n < -1 :=
  1/(a (n + 1)) (-x^(m + 1) (a + b x)^(n + 1) + (m + n + 2) IntegrateTable[x^m (a + b x)^(n + 1), x]);

(* Form 3: Reduces the power of x *)
IntegrateTable[x_^m_ (b_. x_ + a_.)^n_, x_] /; FreeQ[{a, b, m, n}, x] && m + n + 1 =!= 0 && IntegerQ[m] && IntegerQ[n] && m > 0 :=
  1/(b (m + n + 1)) (x^m (a + b x)^(n + 1) - m a IntegrateTable[x^(m - 1) (a + b x)^n, x]);

(* Formula 27 *)
IntegrateTable[1/(a_ + b_. x_), x_] /; FreeQ[{a, b}, x] := 1/b Log[a + b x];

(* Formula 28 *)
IntegrateTable[1/(a_ + b_. x_)^2, x_] /; FreeQ[{a, b}, x] := -1/(b (a + b x));

(* Formula 29 *)
IntegrateTable[1/(a_ + b_. x_)^3, x_] /; FreeQ[{a, b}, x] := -1/(2 b (a + b x)^2);

(* Formula 30 *)
IntegrateTable[x_/(a_ + b_. x_), x_] /; FreeQ[{a, b}, x] := 1/b^2 (a + b x - a Log[a + b x]);
  (* Note: The alternative form x/b - a/b^2 Log[a + b x] is mathematically equivalent *)

(* Formula 31 *)
IntegrateTable[x_/(a_ + b_. x_)^2, x_] /; FreeQ[{a, b}, x] := 1/b^2 (Log[a + b x] + a/(a + b x));

(* Formula 32 *)
IntegrateTable[x_/(a_ + b_. x_)^n_, x_] /; FreeQ[{a, b, n}, x] && n =!= 1 && n =!= 2 := 1/b^2 (-1/((n - 2) (a + b x)^(n - 2)) + a/((n - 1) (a + b x)^(n - 1)));

(* Formula 33 *)
IntegrateTable[x_^2/(a_ + b_. x_), x_] /; FreeQ[{a, b}, x] := 1/b^3 (1/2 (a + b x)^2 - 2 a (a + b x) + a^2 Log[a + b x]);

(* Formula 34 *)
IntegrateTable[x_^2/(a_ + b_. x_)^2, x_] /; FreeQ[{a, b}, x] := 1/b^3 (a + b x - 2 a Log[a + b x] - a^2/(a + b x));

(* Formula 35 *)
IntegrateTable[x_^2/(a_ + b_. x_)^3, x_] /; FreeQ[{a, b}, x] := 1/b^3 (Log[a + b x] + 2 a/(a + b x) - a^2/(2 (a + b x)^2));

(* Formula 36 *)
IntegrateTable[x_^2/(a_ + b_. x_)^n_, x_] /; FreeQ[{a, b, n}, x] && n =!= 1 && n =!= 2 && n =!= 3 := 1/b^3 (-1/((n - 3) (a + b x)^(n - 3)) + 2 a/((n - 2) (a + b x)^(n - 2)) - a^2/((n - 1) (a + b x)^(n - 1)));

(* Formula 37 *)
IntegrateTable[1/(x_ (a_ + b_. x_)), x_] /; FreeQ[{a, b}, x] := 1/a Log[x/(a + b x)];

(* Formula 38 *)
IntegrateTable[1/(x_ (a_ + b_. x_)^2), x_] /; FreeQ[{a, b}, x] := 1/(a (a + b x)) - 1/a^2 Log[(a + b x)/x];

(* Formula 39 *)
IntegrateTable[1/(x_ (a_ + b_. x_)^3), x_] /; FreeQ[{a, b}, x] := 1/a^3 (1/2 ((2 a + b x)/(a + b x))^2 - Log[(a + b x)/x]);

(* Formula 40 *)
IntegrateTable[1/(x_^2 (a_ + b_. x_)), x_] /; FreeQ[{a, b}, x] := -1/(a x) + b/a^2 Log[(a + b x)/x];

(* Formula 41 *)
IntegrateTable[1/(x_^3 (a_ + b_. x_)), x_] /; FreeQ[{a, b}, x] := (2 b x - a)/(2 a^2 x^2) + b^2/a^3 Log[x/(a + b x)];

(* Formula 42 *)
IntegrateTable[1/(x_^2 (a_ + b_. x_)^2), x_] /; FreeQ[{a, b}, x] := -(a + 2 b x)/(a^2 x (a + b x)) + (2 b)/a^3 Log[(a + b x)/x];

(* --- Section 5.4.3: Forms containing c^2 ± x^2 and x^2 - c^2 --- *)

(* Formulas 43-51 (c_^2 ± x_^2, x_^2 - c_^2 and their powers) are omitted:
   with a squared pattern constant `c_^2` they never matched a numeric
   argument, and every one of their integrands is already handled — for
   symbolic and numeric constants alike — by the linear-coefficient forms
   `1/(a_. + b_. x_^2)` (Formulas 60/61), `x_/(a_. + b_. x_^2)` (63),
   `1/(a_ + b_. x_^2)^m_` (67) and `x_/(a_ + b_. x_^2)^m_` (68) below, which
   bind the constant term as a_ and so fire for e.g. 1/(9 + x^2) or
   1/(x^2 - 9)^2 directly. *)

(* --- Section 5.4.4: Forms containing a + bx and c + dx --- *)

(* Formula 52 *)
IntegrateTable[1/((a_ + b_. x_) (c_ + d_. x_)), x_] /; FreeQ[{a, b, c, d}, x] := With[{k = a d - b c, u = a + b x, v = c + d x}, 1/k Log[v/u]];

(* Formula 53 *)
IntegrateTable[x_/((a_ + b_. x_) (c_ + d_. x_)), x_] /; FreeQ[{a, b, c, d}, x] := With[{k = a d - b c, u = a + b x, v = c + d x}, 1/k (a/b Log[u] - c/d Log[v])];

(* Formula 54 *)
IntegrateTable[1/((a_ + b_. x_)^2 (c_ + d_. x_)), x_] /; FreeQ[{a, b, c, d}, x] := With[{k = a d - b c, u = a + b x, v = c + d x}, 1/k (1/u + d/k Log[v/u])];

(* Formula 55 *)
IntegrateTable[x_/((a_ + b_. x_)^2 (c_ + d_. x_)), x_] /; FreeQ[{a, b, c, d}, x] := With[{k = a d - b c, u = a + b x, v = c + d x}, -(a/(b k u)) - c/k^2 Log[v/u]];

(* Formula 56 *)
IntegrateTable[x_^2/((a_ + b_. x_)^2 (c_ + d_. x_)), x_] /; FreeQ[{a, b, c, d}, x] := With[{k = a d - b c, u = a + b x, v = c + d x}, a^2/(b^2 k u) + 1/k^2 (c^2/d Log[v] + (a (k - b c))/b^2 Log[u])];

(* Formula 57 *)
IntegrateTable[1/((a_ + b_. x_)^n_ (c_ + d_. x_)^m_), x_] /; FreeQ[{a, b, c, d, n, m}, x] && m =!= 1 && IntegerQ[m] && m > 1 := With[{k = a d - b c, u = a + b x, v = c + d x}, 1/(k (m - 1)) (-1/(u^(n - 1) v^(m - 1)) - b (m + n - 2) IntegrateTable[1/(u^n v^(m - 1)), x])];
  
(* Formula 58 *)
IntegrateTable[(a_ + b_. x_)/(c_ + d_. x_), x_] /; FreeQ[{a, b, c, d}, x] := With[{k = a d - b c, u = a + b x, v = c + d x}, (b x)/d + k/d^2 Log[v]];


(* Formula 59 *)
IntegrateTable[(a_. + b_. x_)^m_. (c_. + d_. x_)^nn_., x_] /; FreeQ[{a, b, c, d, m, nn}, x] && IntegerQ[nn] && nn < -1 := With[{k = a d - b c, u = a + b x, v = c + d x, n = -nn}, -1/(k (n - 1)) (u^(m + 1)/v^(n - 1) + b (n - m - 2) IntegrateTable[u^m/v^(n - 1), x])];

(* Formula 60 *)
IntegrateTable[(a_. + b_. x_^2)^-1, x_] /; FreeQ[{a, b}, x] && a b > 0 := 1/Sqrt[a b] ArcTan[(x Sqrt[a b])/a];

(* Formula 61 *)
IntegrateTable[(a_. + b_. x_^2)^-1, x_] /; FreeQ[{a, b}, x] && a b < 0 := 1/(2 Sqrt[-a b]) Log[(a + x Sqrt[-a b])/(a - x Sqrt[-a b])];

(* Formula 62 ((a_^2 + b_^2 x_^2)^-1) omitted: squared pattern constants never
   matched a numeric argument, and the integrand is already covered by
   Formula 60, `(a_. + b_. x_^2)^-1` with a b > 0. *)

(* Formula 63 *)
IntegrateTable[x_ (a_. + b_. x_^2)^-1, x_] /; FreeQ[{a, b}, x] := 1/(2 b) Log[a + b x^2];

(* Formula 64 *)
IntegrateTable[x_^2 (a_. + b_. x_^2)^-1, x_] /; FreeQ[{a, b}, x] := x/b - a/b IntegrateTable[(a + b x^2)^-1, x];

(* Formula 65 *)
IntegrateTable[1/(a_ + b_. x_^2)^2, x_] /; FreeQ[{a, b}, x] := x/(2 a (a + b x^2)) + 1/(2 a) IntegrateTable[1/(a + b x^2), x];

(* Formula 66 (1/(a_^2 - b_^2 x_^2)) omitted: squared pattern constants never
   matched a numeric argument, and the integrand is already covered by
   Formula 61, `(a_. + b_. x_^2)^-1` with a b < 0. *)

(* Formula 67: Mapped n+1 to m *)
IntegrateTable[1/(a_ + b_. x_^2)^m_, x_] /; FreeQ[{a, b, m}, x] && m =!= 1 && IntegerQ[m] && m > 1 := x/(2 a (m - 1) (a + b x^2)^(m - 1)) + (2 m - 3)/(2 a (m - 1)) IntegrateTable[1/(a + b x^2)^(m - 1), x];

(* Formula 68: Mapped m+1 to m *)
IntegrateTable[x_/(a_ + b_. x_^2)^m_, x_] /; FreeQ[{a, b, m}, x] && m =!= 1 := -1/(2 b (m - 1) (a + b x^2)^(m - 1));

(* Formula 69: Mapped m+1 to m *)
IntegrateTable[x_^2/(a_ + b_. x_^2)^m_, x_] /; FreeQ[{a, b, m}, x] && m =!= 1 && IntegerQ[m] && m > 1 := -x/(2 b (m - 1) (a + b x^2)^(m - 1)) + 1/(2 b (m - 1)) IntegrateTable[1/(a + b x^2)^(m - 1), x];

(* Formula 70 *)
IntegrateTable[1/(x_ (a_ + b_. x_^2)), x_] /; FreeQ[{a, b}, x] := 1/(2 a) Log[x^2/(a + b x^2)];

(* Formula 71 *)
IntegrateTable[1/(x_^2 (a_ + b_. x_^2)), x_] /; FreeQ[{a, b}, x] := -1/(a x) - b/a IntegrateTable[1/(a + b x^2), x];

(* Formula 72: Mapped n+1 to n *)
IntegrateTable[1/(x_ (a_ + b_. x_^2)^n_), x_] /; FreeQ[{a, b, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(2 a (n - 1) (a + b x^2)^(n - 1)) + 1/a IntegrateTable[1/(x (a + b x^2)^(n - 1)), x];

(* Formula 73: Mapped m+1 to n *)
IntegrateTable[1/(x_^2 (a_ + b_. x_^2)^n_), x_] /; FreeQ[{a, b, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/a IntegrateTable[1/(x^2 (a + b x^2)^(n - 1)), x] - b/a IntegrateTable[1/(a + b x^2)^n, x];


(* Formula 74 *)
IntegrateTable[1/(a_ + b_. x_^3), x_] /; FreeQ[{a, b}, x] := With[{k = (a/b)^(1/3)}, 
    k/(3 a) (1/2 Log[(k + x)^3/(a + b x^3)] + Sqrt[3] ArcTan[(2 x - k)/(k Sqrt[3])])];

(* Formula 75 *)
IntegrateTable[x_/(a_ + b_. x_^3), x_] /; FreeQ[{a, b}, x] := With[{k = (a/b)^(1/3)}, 
    1/(3 b k) (1/2 Log[(a + b x^3)/(k + x)^3] + Sqrt[3] ArcTan[(2 x - k)/(k Sqrt[3])])];

(* Formula 76 *)
IntegrateTable[x_^2/(a_ + b_. x_^3), x_] /; FreeQ[{a, b}, x] := 1/(3 b) Log[a + b x^3];

(* Formula 77: Split into sign conditions *)
IntegrateTable[1/(a_ + b_. x_^4), x_] /; FreeQ[{a, b}, x] && a b > 0 := With[{k = (a/(4 b))^(1/4)}, 
    k/(2 a) (1/2 Log[(x^2 + 2 k x + 2 k^2)/(x^2 - 2 k x + 2 k^2)] + ArcTan[(2 k x)/(2 k^2 - x^2)])];

IntegrateTable[1/(a_ + b_. x_^4), x_] /; FreeQ[{a, b}, x] && a b < 0 := With[{k = (-a/b)^(1/4)}, 
    k/(2 a) (1/2 Log[(x + k)/(x - k)] + ArcTan[x/k])];

(* Formula 78 & 79 *)
IntegrateTable[x_/(a_ + b_. x_^4), x_] /; FreeQ[{a, b}, x] && a b > 0 := With[{k = Sqrt[a/b]}, 1/(2 b k) ArcTan[x^2/k]];

IntegrateTable[x_/(a_ + b_. x_^4), x_] /; FreeQ[{a, b}, x] && a b < 0 := With[{k = Sqrt[-a/b]}, 1/(4 b k) Log[(x^2 - k)/(x^2 + k)]];

(* Formula 80 & 81 *)
IntegrateTable[x_^2/(a_ + b_. x_^4), x_] /; FreeQ[{a, b}, x] && a b > 0 := With[{k = (a/(4 b))^(1/4)}, 
    1/(4 b k) (1/2 Log[(x^2 - 2 k x + 2 k^2)/(x^2 + 2 k x + 2 k^2)] + ArcTan[(2 k x)/(2 k^2 - x^2)])];

IntegrateTable[x_^2/(a_ + b_. x_^4), x_] /; FreeQ[{a, b}, x] && a b < 0 := With[{k = (-a/b)^(1/4)}, 
    1/(4 b k) (Log[(x - k)/(x + k)] + 2 ArcTan[x/k])];

(* Formula 82 *)
IntegrateTable[x_^3/(a_ + b_. x_^4), x_] /; FreeQ[{a, b}, x] := 1/(4 b) Log[a + b x^4];

(* Formula 83 *)
IntegrateTable[1/(x_ (a_ + b_. x_^n_)), x_] /; FreeQ[{a, b, n}, x] && n =!= 0 := 1/(a n) Log[x^n/(a + b x^n)];

(* Formula 84: Mapped m+1 to m *)
IntegrateTable[1/(a_ + b_. x_^n_)^m_, x_] /; FreeQ[{a, b, n, m}, x] && m =!= 1 && IntegerQ[m] && m > 1 := 1/a IntegrateTable[1/(a + b x^n)^(m - 1), x] - b/a IntegrateTable[x^n/(a + b x^n)^m, x];

(* Formula 85: Mapped p+1 to p *)
IntegrateTable[x_^m_/(a_ + b_. x_^n_)^p_, x_] /; FreeQ[{a, b, n, m, p}, x] && p =!= 1 && IntegerQ[p] && p > 1 := 1/b IntegrateTable[x^(m - n)/(a + b x^n)^(p - 1), x] - a/b IntegrateTable[x^(m - n)/(a + b x^n)^p, x];

(* Formula 86: Mapped p+1 to p *)
IntegrateTable[1/(x_^m_ (a_ + b_. x_^n_)^p_), x_] /; FreeQ[{a, b, n, m, p}, x] && p =!= 1 && IntegerQ[p] && p > 1 := 1/a IntegrateTable[1/(x^m (a + b x^n)^(p - 1)), x] - b/a IntegrateTable[1/(x^(m - n) (a + b x^n)^p), x];

(* Formula 87: Form 1 *)
IntegrateTable[x_^m_ (a_ + b_. x_^n_)^p_, x_] /; FreeQ[{a, b, n, m, p}, x] && IntegerQ[m] && IntegerQ[n] && n > 0 && m >= n && n p + m + 1 =!= 0 := 1/(b (n p + m + 1)) (x^(m - n + 1) (a + b x^n)^(p + 1) - a (m - n + 1) IntegrateTable[x^(m - n) (a + b x^n)^p, x]);

(* Formulas 88-106 (c^3 ± x^3, c^4 ± x^4 and their powers).  Squared/cubed
   pattern constants c_^3, c_^4 never matched a numeric argument.  The
   single-power members (88, 91, 94, 96, 101-106) are already handled by the
   linear-coefficient forms above — Formula 74 `1/(a_ + b_. x_^3)`, 75, 76,
   83 `1/(x_ (a_ + b_. x_^n_))`, and 77-82 for x^4 — which bind the constant as
   a_ and so fire for e.g. 1/(8 + x^3), x/(16 + x^4).  The higher-power
   reductions below are kept, with the cubed constant bound linearly as c_
   (c^3 -> c, c^6 -> c^2); their recursion bottoms out in the linear forms. *)

(* Formula 89 *)
IntegrateTable[1/(c_ + x_^3)^2, x_] /; FreeQ[c, x] := x/(3 c (c + x^3)) + 2/(3 c) IntegrateTable[1/(c + x^3), x];
IntegrateTable[1/(c_ - x_^3)^2, x_] /; FreeQ[c, x] := x/(3 c (c - x^3)) + 2/(3 c) IntegrateTable[1/(c - x^3), x];

(* Formula 90 *)
IntegrateTable[1/(c_ + x_^3)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(3 (n - 1) c) (x/(c + x^3)^(n - 1) + (3 n - 4) IntegrateTable[1/(c + x^3)^(n - 1), x]);
IntegrateTable[1/(c_ - x_^3)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(3 (n - 1) c) (x/(c - x^3)^(n - 1) + (3 n - 4) IntegrateTable[1/(c - x^3)^(n - 1), x]);

(* Formula 92 *)
IntegrateTable[x_/(c_ + x_^3)^2, x_] /; FreeQ[c, x] := x^2/(3 c (c + x^3)) + 1/(3 c) IntegrateTable[x/(c + x^3), x];
IntegrateTable[x_/(c_ - x_^3)^2, x_] /; FreeQ[c, x] := x^2/(3 c (c - x^3)) + 1/(3 c) IntegrateTable[x/(c - x^3), x];

(* Formula 93 *)
IntegrateTable[x_/(c_ + x_^3)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(3 (n - 1) c) (x^2/(c + x^3)^(n - 1) + (3 n - 5) IntegrateTable[x/(c + x^3)^(n - 1), x]);
IntegrateTable[x_/(c_ - x_^3)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(3 (n - 1) c) (x^2/(c - x^3)^(n - 1) + (3 n - 5) IntegrateTable[x/(c - x^3)^(n - 1), x]);

(* Formula 95.  (An explicit n = 2 form is given alongside the general n_ rule:
   a pattern exponent `(c_ + x_^3)^n_` in a denominator does not match a
   literal negative power, whereas the explicit `^2` canonicalises and does.) *)
IntegrateTable[x_^2/(c_ + x_^3)^2, x_] /; FreeQ[c, x] := -1/(3 (c + x^3));
IntegrateTable[x_^2/(c_ - x_^3)^2, x_] /; FreeQ[c, x] := 1/(3 (c - x^3));
IntegrateTable[x_^2/(c_ + x_^3)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 := -1/(3 (n - 1) (c + x^3)^(n - 1));
IntegrateTable[x_^2/(c_ - x_^3)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 := 1/(3 (n - 1) (c - x^3)^(n - 1));

(* Formula 97 *)
IntegrateTable[1/(x_ (c_ + x_^3)^2), x_] /; FreeQ[c, x] := 1/(3 c (c + x^3)) + 1/(3 c^2) Log[x^3/(c + x^3)];
IntegrateTable[1/(x_ (c_ - x_^3)^2), x_] /; FreeQ[c, x] := 1/(3 c (c - x^3)) + 1/(3 c^2) Log[x^3/(c - x^3)];

(* Formula 98 *)
IntegrateTable[1/(x_ (c_ + x_^3)^n_), x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(3 (n - 1) c (c + x^3)^(n - 1)) + 1/c IntegrateTable[1/(x (c + x^3)^(n - 1)), x];
IntegrateTable[1/(x_ (c_ - x_^3)^n_), x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(3 (n - 1) c (c - x^3)^(n - 1)) + 1/c IntegrateTable[1/(x (c - x^3)^(n - 1)), x];

(* Formula 99 *)
IntegrateTable[1/(x_^2 (c_ + x_^3)), x_] /; FreeQ[c, x] := -1/(c x) - 1/c IntegrateTable[x/(c + x^3), x];
IntegrateTable[1/(x_^2 (c_ - x_^3)), x_] /; FreeQ[c, x] := -1/(c x) + 1/c IntegrateTable[x/(c - x^3), x];

(* Formula 100 *)
IntegrateTable[1/(x_^2 (c_ + x_^3)^n_), x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/c IntegrateTable[1/(x^2 (c + x^3)^(n - 1)), x] - 1/c IntegrateTable[x/(c + x^3)^n, x];
IntegrateTable[1/(x_^2 (c_ - x_^3)^n_), x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/c IntegrateTable[1/(x^2 (c - x^3)^(n - 1)), x] + 1/c IntegrateTable[x/(c - x^3)^n, x];

(* Formula 107 (x_^3/(c_^4 ± x_^4)) omitted: the squared pattern constant never
   matched a numeric argument, and the integrand is already covered by
   Formula 82, `x_^3/(a_ + b_. x_^4)`, which binds the constant as a_. *)

(* Formula 108 *)
IntegrateTable[1/(a_ + b_. x_ + c_. x_^2), x_] /; FreeQ[{a, b, c}, x] && 4 a c - b^2 > 0 := With[{q = 4 a c - b^2}, 2/Sqrt[q] ArcTan[(2 c x + b)/Sqrt[q]]];

IntegrateTable[1/(a_ + b_. x_ + c_. x_^2), x_] /; FreeQ[{a, b, c}, x] && 4 a c - b^2 < 0 := With[{q = 4 a c - b^2}, 1/Sqrt[-q] Log[(2 c x + b - Sqrt[-q])/(2 c x + b + Sqrt[-q])]];

(* Formula 109 *)
IntegrateTable[1/(a_ + b_. x_ + c_. x_^2)^2, x_] /; FreeQ[{a, b, c}, x] := With[{X = a + b x + c x^2, q = 4 a c - b^2}, 
    (2 c x + b)/(q X) + (2 c)/q IntegrateTable[1/X, x]];

(* Formula 110 *)
IntegrateTable[1/(a_ + b_. x_ + c_. x_^2)^3, x_] /; FreeQ[{a, b, c}, x] := With[{X = a + b x + c x^2, q = 4 a c - b^2}, 
    (2 c x + b)/q (1/(2 X^2) + (3 c)/(q X)) + (6 c^2)/q^2 IntegrateTable[1/X, x]];

(* Formula 111: Mapped n+1 to n *)
IntegrateTable[1/(a_ + b_. x_ + c_. x_^2)^n_, x_] /; FreeQ[{a, b, c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := With[{X = a + b x + c x^2, q = 4 a c - b^2},
    (2 c x + b)/((n - 1) q X^(n - 1)) + (2 (2 n - 3) c)/(q (n - 1)) IntegrateTable[1/X^(n - 1), x]];

(* Formula 112 *)
IntegrateTable[x_/(a_ + b_. x_ + c_. x_^2), x_] /; FreeQ[{a, b, c}, x] := With[{X = a + b x + c x^2}, 
    1/(2 c) Log[X] - b/(2 c) IntegrateTable[1/X, x]];

(* Formula 113 *)
IntegrateTable[x_/(a_ + b_. x_ + c_. x_^2)^2, x_] /; FreeQ[{a, b, c}, x] := With[{X = a + b x + c x^2, q = 4 a c - b^2}, 
    -(b x + 2 a)/(q X) - b/q IntegrateTable[1/X, x]];

(* Formula 114: Mapped n+1 to n *)
IntegrateTable[x_/(a_ + b_. x_ + c_. x_^2)^n_, x_] /; FreeQ[{a, b, c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := With[{X = a + b x + c x^2, q = 4 a c - b^2},
    -(2 a + b x)/((n - 1) q X^(n - 1)) - (b (2 n - 3))/((n - 1) q) IntegrateTable[1/X^(n - 1), x]];

(* Formula 115 *)
IntegrateTable[x_^2/(a_ + b_. x_ + c_. x_^2), x_] /; FreeQ[{a, b, c}, x] := With[{X = a + b x + c x^2}, 
    x/c - b/(2 c^2) Log[X] + (b^2 - 2 a c)/(2 c^2) IntegrateTable[1/X, x]];

(* Formula 116 *)
IntegrateTable[x_^2/(a_ + b_. x_ + c_. x_^2)^2, x_] /; FreeQ[{a, b, c}, x] := With[{X = a + b x + c x^2, q = 4 a c - b^2}, 
    ((b^2 - 2 a c) x + a b)/(c q X) + (2 a)/q IntegrateTable[1/X, x]];

(* Formula 117: Mapped n+1 to n *)
IntegrateTable[x_^m_/(a_ + b_. x_ + c_. x_^2)^n_, x_] /; FreeQ[{a, b, c, m, n}, x] && n =!= 1 && IntegerQ[m] && m > 1 := With[{X = a + b x + c x^2},
    -x^(m - 1)/((2 n - m - 1) c X^(n - 1)) - (n - m)*b/((2 n - m - 1) c) IntegrateTable[x^(m - 1)/X^n, x] + (m - 1)*a/((2 n - m - 1) c) IntegrateTable[x^(m - 2)/X^n, x]];

(* Formula 118 *)
IntegrateTable[1/(x_ (a_ + b_. x_ + c_. x_^2)), x_] /; FreeQ[{a, b, c}, x] := With[{X = a + b x + c x^2}, 1/(2 a) Log[x^2/X] - b/(2 a) IntegrateTable[1/X, x]];

(* Formula 119 *)
IntegrateTable[1/(x_^2 (a_ + b_. x_ + c_. x_^2)), x_] /; FreeQ[{a, b, c}, x] := With[{X = a + b x + c x^2}, 
    b/(2 a^2) Log[X/x^2] - 1/(a x) + (b^2/(2 a^2) - c/a) IntegrateTable[1/X, x]];

(* Formula 120: Mapped n to n-1 in denominators for smooth matching *)
IntegrateTable[1/(x_ (a_ + b_. x_ + c_. x_^2)^n_), x_] /; FreeQ[{a, b, c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := With[{X = a + b x + c x^2},
    1/(2 a (n - 1) X^(n - 1)) - b/(2 a) IntegrateTable[1/X^n, x] + 1/a IntegrateTable[1/(x X^(n - 1)), x]];

(* Formula 121: Mapped n+1 to n for smoother general matching *)
IntegrateTable[1/(x_^m_ (a_ + b_. x_ + c_. x_^2)^n_), x_] /; FreeQ[{a, b, c, m, n}, x] && m =!= 1 && IntegerQ[m] && m > 1 := With[{X = a + b x + c x^2},
    -1/((m - 1) a x^(m - 1) X^(n - 1)) - ((n + m - 2) b)/((m - 1) a) IntegrateTable[1/(x^(m - 1) X^n), x] - ((2 n + m - 3) c)/((m - 1) a) IntegrateTable[1/(x^(m - 2) X^n), x]];

(* Formula 122 *)
IntegrateTable[Sqrt[a_ + b_. x_], x_] /; FreeQ[{a, b}, x] := 2/(3 b) (a + b x)^(3/2);

(* Formula 123 *)
IntegrateTable[x_ Sqrt[a_ + b_. x_], x_] /; FreeQ[{a, b}, x] := (-2 (2 a - 3 b x))/(15 b^2) (a + b x)^(3/2);

(* Formula 124 *)
IntegrateTable[x_^2 Sqrt[a_ + b_. x_], x_] /; FreeQ[{a, b}, x] := (2 (8 a^2 - 12 a b x + 15 b^2 x^2))/(105 b^3) (a + b x)^(3/2);

(* Formula 125.  A downward recursion m -> m - 1 that terminates at the m = 0
 * base case (Formula 122).  It is valid only for positive integer m: for a
 * non-integer m (e.g. Sqrt[x] Sqrt[a + b x], where m = 1/2) it would march
 * m -> -1/2 -> -3/2 -> ... without ever reaching the base case, recursing until
 * $RecursionLimit and emitting a spurious 1/0 as it crosses m = -3/2 (where the
 * 2 m + 3 denominator vanishes).  Guarded like its siblings (Formulas 127, 134). *)
IntegrateTable[x_^m_ Sqrt[a_ + b_. x_], x_] /; FreeQ[{a, b, m}, x] && IntegerQ[m] && m > 0 := 2/(b (2 m + 3)) (x^m (a + b x)^(3/2) - m a IntegrateTable[x^(m - 1) Sqrt[a + b x], x]);

(* Formula 126 *)
IntegrateTable[Sqrt[a_ + b_. x_]/x_, x_] /; FreeQ[{a, b}, x] := 2 Sqrt[a + b x] + a IntegrateTable[1/(x Sqrt[a + b x]), x];

(* Formula 127 *)
IntegrateTable[Sqrt[a_ + b_. x_]/x_^m_, x_] /; FreeQ[{a, b, m}, x] && m =!= 1 && IntegerQ[m] && m > 1 := -1/((m - 1) a) ((a + b x)^(3/2)/x^(m - 1) + ((2 m - 5) b)/2 IntegrateTable[Sqrt[a + b x]/x^(m - 1), x]);

(* Formula 128 *)
IntegrateTable[1/Sqrt[a_ + b_. x_], x_] /; FreeQ[{a, b}, x] := 2/b Sqrt[a + b x];

(* Formula 129 *)
IntegrateTable[x_/Sqrt[a_ + b_. x_], x_] /; FreeQ[{a, b}, x] := (-2 (2 a - b x))/(3 b^2) Sqrt[a + b x];

(* Formula 130 *)
IntegrateTable[x_^2/Sqrt[a_ + b_. x_], x_] /; FreeQ[{a, b}, x] := (2 (8 a^2 - 4 a b x + 3 b^2 x^2))/(15 b^3) Sqrt[a + b x];

(* Formula 131.  As with Formula 125, this m -> m - 1 recursion terminates at the
 * m = 0 base case (Formula 128) only for positive integer m; a non-integer m
 * would recurse unboundedly and hit a 1/0 at m = -1/2 (vanishing 2 m + 1). *)
IntegrateTable[x_^m_/Sqrt[a_ + b_. x_], x_] /; FreeQ[{a, b, m}, x] && IntegerQ[m] && m > 0 := 2/((2 m + 1) b) (x^m Sqrt[a + b x] - m a IntegrateTable[x^(m - 1)/Sqrt[a + b x], x]);

(* Formula 132 *)
IntegrateTable[1/(x_ Sqrt[a_ + b_. x_]), x_] /; FreeQ[{a, b}, x] && a < 0 := 2/Sqrt[-a] ArcTan[Sqrt[(a + b x)/-a]];

IntegrateTable[1/(x_ Sqrt[a_ + b_. x_]), x_] /; FreeQ[{a, b}, x] && a > 0 := 1/Sqrt[a] Log[(Sqrt[a + b x] - Sqrt[a])/(Sqrt[a + b x] + Sqrt[a])];

(* Formula 133 *)
IntegrateTable[1/(x_^2 Sqrt[a_ + b_. x_]), x_] /; FreeQ[{a, b}, x] := -Sqrt[a + b x]/(a x) - b/(2 a) IntegrateTable[1/(x Sqrt[a + b x]), x];

(* Formula 134 *)
IntegrateTable[1/(x_^n_ Sqrt[a_ + b_. x_]), x_] /; FreeQ[{a, b, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := -Sqrt[a + b x]/((n - 1) a x^(n - 1)) - ((2 n - 3) b)/(2 (n - 1) a) IntegrateTable[1/(x^(n - 1) Sqrt[a + b x]), x];

(* Formula 135 & 136: Converted textbook ±n/2 to general rational power m_ *)
IntegrateTable[(a_ + b_. x_)^m_, x_] /; FreeQ[{a, b, m}, x] && m =!= -1 := (a + b x)^(m + 1)/(b (m + 1));

IntegrateTable[x_ (a_ + b_. x_)^m_, x_] /; FreeQ[{a, b, m}, x] && m =!= -1 && m =!= -2 := 1/b^2 ((a + b x)^(m + 2)/(m + 2) - a (a + b x)^(m + 1)/(m + 1));

(* Formula 137 & 138: Using general rational power m_ *)
IntegrateTable[1/(x_ (a_ + b_. x_)^m_), x_] /; FreeQ[{a, b, m}, x] && IntegerQ[m] && m > 1 := 1/a IntegrateTable[1/(x (a + b x)^(m - 1)), x] - b/a IntegrateTable[1/(a + b x)^m, x];

IntegrateTable[(a_ + b_. x_)^m_/x_, x_] /; FreeQ[{a, b, m}, x] && IntegerQ[m] && m > 0 := b IntegrateTable[(a + b x)^(m - 1), x] + a IntegrateTable[(a + b x)^(m - 1)/x, x];

(* Formula 139 *)
IntegrateTable[1/Sqrt[(a_ + b_. x_) (c_ + d_. x_)], x_] /; FreeQ[{a, b, c, d}, x] && b d > 0 := 2/Sqrt[b d] ArcTanh[Sqrt[b d (a + b x) (c + d x)]/(b (c + d x))];

IntegrateTable[1/Sqrt[(a_ + b_. x_) (c_ + d_. x_)], x_] /; FreeQ[{a, b, c, d}, x] && b d < 0 := 2/Sqrt[-b d] ArcTan[Sqrt[-b d (a + b x) (c + d x)]/(b (c + d x))];

(* Formula 140 *)
IntegrateTable[Sqrt[(a_ + b_. x_) (c_ + d_. x_)], x_] /; FreeQ[{a, b, c, d}, x] := With[{u = a + b x, v = c + d x, k = a d - b c}, 
    (k + 2 b v)/(4 b d) Sqrt[u v] - k^2/(8 b d) IntegrateTable[1/Sqrt[u v], x]];

(* Formula 141 *)
(* Integral[1/((a+b x) Sqrt[(a+b x)(c+d x)]), x] = (2/k) Sqrt[(c+d x)/(a+b x)],
 * where k = a d - b c.  Verified by D[(2/k) Sqrt[v/u], x] = 1/(u Sqrt[u v]).
 * The previous form `(2/(k Sqrt[k d])) ArcTanh[Sqrt[k d u v]/(d u)]` was
 * off by a factor of (a+b x): differentiating it yields 1/Sqrt[u v] rather
 * than 1/(u Sqrt[u v]).  Numerical and symbolic verification:
 *     D[(2/k) Sqrt[v/u], x] - 1/(u Sqrt[u v])  ==  0  (Simplify) *)
IntegrateTable[1/((a_ + b_. x_) Sqrt[(a_ + b_. x_) (c_ + d_. x_)]), x_] /; FreeQ[{a, b, c, d}, x] && (a d - b c) d > 0 := With[{u = a + b x, v = c + d x, k = a d - b c},
    (2/k) Sqrt[v/u]];

(* Formula 142 *)
IntegrateTable[x_/Sqrt[(a_ + b_. x_) (c_ + d_. x_)], x_] /; FreeQ[{a, b, c, d}, x] := With[{u = a + b x, v = c + d x}, 
    Sqrt[u v]/(b d) - (a d + b c)/(2 b d) IntegrateTable[1/Sqrt[u v], x]];

(* Formula 143 *)
IntegrateTable[1/((c_ + d_. x_) Sqrt[(a_ + b_. x_) (c_ + d_. x_)]), x_] /; FreeQ[{a, b, c, d}, x] := With[{u = a + b x, v = c + d x, k = a d - b c}, 
    -(2 Sqrt[u v])/(k v)];

(* Formula 144 *)
IntegrateTable[(c_ + d_. x_)/Sqrt[(a_ + b_. x_) (c_ + d_. x_)], x_] /; FreeQ[{a, b, c, d}, x] := With[{u = a + b x, v = c + d x, k = a d - b c}, 
    Sqrt[u v]/b - k/(2 b) IntegrateTable[1/Sqrt[u v], x]];

(* Formula 145 *)
IntegrateTable[Sqrt[(c_ + d_. x_)/(a_ + b_. x_)], x_] /; FreeQ[{a, b, c, d}, x] := With[{u = a + b x, v = c + d x}, 
    v/Abs[v] IntegrateTable[v/Sqrt[u v], x]];

(* Formula 146 *)
IntegrateTable[(c_ + d_. x_)^m_ Sqrt[a_ + b_. x_], x_] /; FreeQ[{a, b, c, d, m}, x] := With[{u = a + b x, v = c + d x, k = a d - b c}, 
    1/((2 m + 3) d) (2 v^(m + 1) Sqrt[u] + k IntegrateTable[v^m/Sqrt[u], x])];

(* Formula 147 *)
IntegrateTable[1/((c_ + d_. x_)^m_ Sqrt[a_ + b_. x_]), x_] /; FreeQ[{a, b, c, d, m}, x] && m =!= 1 && IntegerQ[m] && m > 1 := With[{u = a + b x, v = c + d x, k = a d - b c},
    -1/((m - 1) k) (Sqrt[u]/v^(m - 1) + (m - 3/2) b IntegrateTable[1/(v^(m - 1) Sqrt[u]), x])];

(* Formula 148 *)
IntegrateTable[(c_ + d_. x_)^m_/Sqrt[a_ + b_. x_], x_] /; FreeQ[{a, b, c, d, m}, x] && IntegerQ[m] && m > 0 := With[{u = a + b x, v = c + d x, k = a d - b c},
    2/(b (2 m + 1)) (v^m Sqrt[u] - m k IntegrateTable[v^(m - 1)/Sqrt[u], x])];

(* Formulas 149-190.  The squared constant a_^2 (which the matcher cannot
   bind against a numeric argument) is replaced by a single a_ bound to the
   constant, with the linear part re-expressed via Sqrt[a].  Both signs are
   matched by the `x_^2 + a_` form: when the RHS is even in the original
   constant a single unguarded rule reproduces the `x^2 - a^2` case for
   a < 0 automatically; when odd (Abs[a], ArcSec, 1/a) the plus and minus
   cases are split by `Not[TrueQ[a < 0]]` / `TrueQ[a < 0]`. *)

(* Formula 149 *)
IntegrateTable[Sqrt[x_^2 + a_], x_] /; FreeQ[a, x] := 1/2 (x Sqrt[x^2 + a] + a Log[x + Sqrt[x^2 + a]]);

(* Formula 150 *)
IntegrateTable[1/Sqrt[x_^2 + a_], x_] /; FreeQ[a, x] := Log[x + Sqrt[x^2 + a]];

(* Formula 152 *)
IntegrateTable[1/(x_ Sqrt[x_^2 + a_]), x_] /; FreeQ[a, x] && Not[TrueQ[a < 0]] := -1/Sqrt[a] Log[(Sqrt[a] + Sqrt[x^2 + a])/x];
(* Formula 151 *)
IntegrateTable[1/(x_ Sqrt[x_^2 + a_]), x_] /; FreeQ[a, x] && TrueQ[a < 0] := 1/Sqrt[-a] ArcSec[x/Sqrt[-a]];

(* Formula 153 *)
IntegrateTable[Sqrt[x_^2 + a_]/x_, x_] /; FreeQ[a, x] && Not[TrueQ[a < 0]] := Sqrt[x^2 + a] - Sqrt[a] Log[(Sqrt[a] + Sqrt[x^2 + a])/x];
(* Formula 154 *)
IntegrateTable[Sqrt[x_^2 + a_]/x_, x_] /; FreeQ[a, x] && TrueQ[a < 0] := Sqrt[x^2 + a] - Sqrt[-a] ArcSec[x/Sqrt[-a]];

(* Formula 155 *)
IntegrateTable[x_/Sqrt[x_^2 + a_], x_] /; FreeQ[a, x] := Sqrt[x^2 + a];

(* Formula 156 *)
IntegrateTable[x_ Sqrt[x_^2 + a_], x_] /; FreeQ[a, x] := 1/3 (x^2 + a)^(3/2);

(* Formula 157 *)
IntegrateTable[(x_^2 + a_)^(3/2), x_] /; FreeQ[a, x] := 1/4 (x (x^2 + a)^(3/2) + (3 a x)/2 Sqrt[x^2 + a] + (3 a^2)/2 Log[x + Sqrt[x^2 + a]]);

(* Formula 158 *)
IntegrateTable[1/(x_^2 + a_)^(3/2), x_] /; FreeQ[a, x] := x/(a Sqrt[x^2 + a]);

(* Formula 159 *)
IntegrateTable[x_/(x_^2 + a_)^(3/2), x_] /; FreeQ[a, x] := -1/Sqrt[x^2 + a];

(* Formula 160 *)
IntegrateTable[x_ (x_^2 + a_)^(3/2), x_] /; FreeQ[a, x] := 1/5 (x^2 + a)^(5/2);

(* Formula 161 *)
IntegrateTable[x_^2 Sqrt[x_^2 + a_], x_] /; FreeQ[a, x] := x/4 (x^2 + a)^(3/2) - (a x)/8 Sqrt[x^2 + a] - a^2/8 Log[x + Sqrt[x^2 + a]];

(* Formula 162 & 163 *)
IntegrateTable[x_^3 Sqrt[x_^2 + a_], x_] /; FreeQ[a, x] := 1/15 (3 x^2 - 2 a) (x^2 + a)^(3/2);

(* Formula 164 *)
IntegrateTable[x_^2/Sqrt[x_^2 + a_], x_] /; FreeQ[a, x] := x/2 Sqrt[x^2 + a] - a/2 Log[x + Sqrt[x^2 + a]];

(* Formula 165 *)
IntegrateTable[x_^3/Sqrt[x_^2 + a_], x_] /; FreeQ[a, x] := 1/3 (x^2 + a)^(3/2) - a Sqrt[x^2 + a];

(* Formula 166 *)
IntegrateTable[1/(x_^2 Sqrt[x_^2 + a_]), x_] /; FreeQ[a, x] := -Sqrt[x^2 + a]/(a x);

(* Formula 167 *)
IntegrateTable[1/(x_^3 Sqrt[x_^2 + a_]), x_] /; FreeQ[a, x] && Not[TrueQ[a < 0]] := -Sqrt[x^2 + a]/(2 a x^2) + 1/(2 a Sqrt[a]) Log[(Sqrt[a] + Sqrt[x^2 + a])/x];
(* Formula 168 *)
IntegrateTable[1/(x_^3 Sqrt[x_^2 + a_]), x_] /; FreeQ[a, x] && TrueQ[a < 0] := -Sqrt[x^2 + a]/(2 a x^2) + 1/(2 (-a) Sqrt[-a]) ArcSec[x/Sqrt[-a]];

(* Formula 169 *)
IntegrateTable[x_^2 (x_^2 + a_)^(3/2), x_] /; FreeQ[a, x] := x/6 (x^2 + a)^(5/2) - (a x)/24 (x^2 + a)^(3/2) - (a^2 x)/16 Sqrt[x^2 + a] - a^3/16 Log[x + Sqrt[x^2 + a]];

(* Formula 170 *)
IntegrateTable[x_^3 (x_^2 + a_)^(3/2), x_] /; FreeQ[a, x] := 1/7 (x^2 + a)^(7/2) - a/5 (x^2 + a)^(5/2);

(* Formula 171 *)
IntegrateTable[Sqrt[x_^2 + a_]/x_^2, x_] /; FreeQ[a, x] := -Sqrt[x^2 + a]/x + Log[x + Sqrt[x^2 + a]];

(* Formula 172 *)
IntegrateTable[Sqrt[x_^2 + a_]/x_^3, x_] /; FreeQ[a, x] && Not[TrueQ[a < 0]] := -Sqrt[x^2 + a]/(2 x^2) - 1/(2 Sqrt[a]) Log[(Sqrt[a] + Sqrt[x^2 + a])/x];
(* Formula 173 *)
IntegrateTable[Sqrt[x_^2 + a_]/x_^3, x_] /; FreeQ[a, x] && TrueQ[a < 0] := -Sqrt[x^2 + a]/(2 x^2) + 1/(2 Sqrt[-a]) ArcSec[x/Sqrt[-a]];

(* Formula 174 *)
IntegrateTable[Sqrt[x_^2 + a_]/x_^4, x_] /; FreeQ[a, x] := -(x^2 + a)^(3/2)/(3 a x^3);

(* Formula 175 *)
IntegrateTable[x_^2/(x_^2 + a_)^(3/2), x_] /; FreeQ[a, x] := -x/Sqrt[x^2 + a] + Log[x + Sqrt[x^2 + a]];

(* Formula 176 *)
IntegrateTable[x_^3/(x_^2 + a_)^(3/2), x_] /; FreeQ[a, x] := Sqrt[x^2 + a] + a/Sqrt[x^2 + a];

(* Formula 177 *)
IntegrateTable[1/(x_ (x_^2 + a_)^(3/2)), x_] /; FreeQ[a, x] && Not[TrueQ[a < 0]] := 1/(a Sqrt[x^2 + a]) - 1/(a Sqrt[a]) Log[(Sqrt[a] + Sqrt[x^2 + a])/x];
(* Formula 178 *)
IntegrateTable[1/(x_ (x_^2 + a_)^(3/2)), x_] /; FreeQ[a, x] && TrueQ[a < 0] := 1/(a Sqrt[x^2 + a]) - 1/((-a) Sqrt[-a]) ArcSec[x/Sqrt[-a]];

(* Formula 179 *)
IntegrateTable[1/(x_^2 (x_^2 + a_)^(3/2)), x_] /; FreeQ[a, x] := -1/a^2 (Sqrt[x^2 + a]/x + x/Sqrt[x^2 + a]);

(* Formula 180 *)
IntegrateTable[1/(x_^3 (x_^2 + a_)^(3/2)), x_] /; FreeQ[a, x] && Not[TrueQ[a < 0]] := -1/(2 a x^2 Sqrt[x^2 + a]) - 3/(2 a^2 Sqrt[x^2 + a]) + 3/(2 a^2 Sqrt[a]) Log[(Sqrt[a] + Sqrt[x^2 + a])/x];
(* Formula 181.  (The CRC reference prints the first term as
   Sqrt[x^2-a^2]/(2 a^2 x^2); that is a transcription error — mirroring the
   x^2+a^2 sibling 180 and re-deriving gives 1/(2 a^2 x^2 Sqrt[...]).) *)
IntegrateTable[1/(x_^3 (x_^2 + a_)^(3/2)), x_] /; FreeQ[a, x] && TrueQ[a < 0] :=
  -1/(2 a x^2 Sqrt[x^2 + a]) - 3/(2 a^2 Sqrt[x^2 + a]) - 3/(2 a^2 Sqrt[-a]) ArcSec[x/Sqrt[-a]];

(* Formula 182: Recursive reductions for positive powers *)
IntegrateTable[x_^m_/Sqrt[x_^2 + a_], x_] /; FreeQ[{a, m}, x] && m =!= 0 && IntegerQ[m] && m > 1 := x^(m - 1)/m Sqrt[x^2 + a] - ((m - 1) a)/m IntegrateTable[x^(m - 2)/Sqrt[x^2 + a], x];

(* Formula 185: Recursive reductions for negative powers *)
IntegrateTable[1/(x_^m_ Sqrt[x_^2 + a_]), x_] /; FreeQ[{a, m}, x] && m =!= 1 && IntegerQ[m] && m > 2 := -Sqrt[x^2 + a]/((m - 1) a x^(m - 1)) - (m - 2)/((m - 1) a) IntegrateTable[1/(x^(m - 2) Sqrt[x^2 + a]), x];

(* Formula 189 & 190.  Here a_ is bound linearly from (x -+ a_), but the
   matcher still cannot verify the `a_^2` inside the radical, so the radical
   constant is matched separately as b_ (= -a^2) and the plus form is used. *)
IntegrateTable[1/((x_ - a_) Sqrt[x_^2 + b_]), x_] /; FreeQ[{a, b}, x] && b === -a^2 := -Sqrt[x^2 + b]/(a (x - a));
IntegrateTable[1/((x_ + a_) Sqrt[x_^2 + b_]), x_] /; FreeQ[{a, b}, x] && b === -a^2 := Sqrt[x^2 + b]/(a (x + a));

(* Formulas 191-220.  The constant a_^2 is always positive here (a^2 - x^2),
   so it is bound as a single a_ (= a^2 > 0) and the linear part recovered as
   Sqrt[a]; `a_ - x_^2` binds cleanly against a numeric constant and no sign
   guard is needed.  a^2 -> a, a^4 -> a^2, a^6 -> a^3, Abs[a] -> Sqrt[a]. *)

(* Formula 191 *)
IntegrateTable[Sqrt[a_ - x_^2], x_] /; FreeQ[a, x] :=
  1/2 (x Sqrt[a - x^2] + a ArcSin[x/Sqrt[a]]);

(* Formula 192 *)
IntegrateTable[1/Sqrt[a_ - x_^2], x_] /; FreeQ[a, x] := ArcSin[x/Sqrt[a]];

(* Formula 193 *)
IntegrateTable[1/(x_ Sqrt[a_ - x_^2]), x_] /; FreeQ[a, x] := -1/Sqrt[a] Log[(Sqrt[a] + Sqrt[a - x^2])/x];

(* Formula 194 *)
IntegrateTable[Sqrt[a_ - x_^2]/x_, x_] /; FreeQ[a, x] := Sqrt[a - x^2] - Sqrt[a] Log[(Sqrt[a] + Sqrt[a - x^2])/x];

(* Formula 195 & 196 *)
IntegrateTable[x_/Sqrt[a_ - x_^2], x_] /; FreeQ[a, x] := -Sqrt[a - x^2];
IntegrateTable[x_ Sqrt[a_ - x_^2], x_] /; FreeQ[a, x] := -1/3 (a - x^2)^(3/2);

(* Formula 197 *)
IntegrateTable[(a_ - x_^2)^(3/2), x_] /; FreeQ[a, x] :=
  1/4 (x (a - x^2)^(3/2) + (3 a x)/2 Sqrt[a - x^2] + (3 a^2)/2 ArcSin[x/Sqrt[a]]);

(* Formula 198 & 199 *)
IntegrateTable[1/(a_ - x_^2)^(3/2), x_] /; FreeQ[a, x] := x/(a Sqrt[a - x^2]);
IntegrateTable[x_/(a_ - x_^2)^(3/2), x_] /; FreeQ[a, x] := 1/Sqrt[a - x^2];

(* Formula 200 *)
IntegrateTable[x_ (a_ - x_^2)^(3/2), x_] /; FreeQ[a, x] := -1/5 (a - x^2)^(5/2);

(* Formula 201 *)
IntegrateTable[x_^2 Sqrt[a_ - x_^2], x_] /; FreeQ[a, x] :=
  -x/4 (a - x^2)^(3/2) + a/8 (x Sqrt[a - x^2] + a ArcSin[x/Sqrt[a]]);

(* Formula 202 *)
IntegrateTable[x_^3 Sqrt[a_ - x_^2], x_] /; FreeQ[a, x] := (-1/5 x^2 - 2/15 a) (a - x^2)^(3/2);

(* Formula 203 *)
IntegrateTable[x_^2 (a_ - x_^2)^(3/2), x_] /; FreeQ[a, x] :=
  -1/6 x (a - x^2)^(5/2) + (a x)/24 (a - x^2)^(3/2) + (a^2 x)/16 Sqrt[a - x^2] + a^3/16 ArcSin[x/Sqrt[a]];

(* Formula 204 *)
IntegrateTable[x_^3 (a_ - x_^2)^(3/2), x_] /; FreeQ[a, x] := 1/7 (a - x^2)^(7/2) - a/5 (a - x^2)^(5/2);

(* Formula 205 *)
IntegrateTable[x_^2/Sqrt[a_ - x_^2], x_] /; FreeQ[a, x] :=
  -x/2 Sqrt[a - x^2] + a/2 ArcSin[x/Sqrt[a]];

(* Formula 206 *)
IntegrateTable[1/(x_^2 Sqrt[a_ - x_^2]), x_] /; FreeQ[a, x] := -Sqrt[a - x^2]/(a x);

(* Formula 207 *)
IntegrateTable[Sqrt[a_ - x_^2]/x_^2, x_] /; FreeQ[a, x] := -Sqrt[a - x^2]/x - ArcSin[x/Sqrt[a]];

(* Formula 208 *)
IntegrateTable[Sqrt[a_ - x_^2]/x_^3, x_] /; FreeQ[a, x] :=
  -Sqrt[a - x^2]/(2 x^2) + 1/(2 Sqrt[a]) Log[(Sqrt[a] + Sqrt[a - x^2])/x];

(* Formula 209 *)
IntegrateTable[Sqrt[a_ - x_^2]/x_^4, x_] /; FreeQ[a, x] := -(a - x^2)^(3/2)/(3 a x^3);

(* Formula 210 *)
IntegrateTable[x_^2/(a_ - x_^2)^(3/2), x_] /; FreeQ[a, x] := x/Sqrt[a - x^2] - ArcSin[x/Sqrt[a]];

(* Formula 211 *)
IntegrateTable[x_^3/Sqrt[a_ - x_^2], x_] /; FreeQ[a, x] := -2/3 (a - x^2)^(3/2) - x^2 Sqrt[a - x^2];

(* Formula 212 *)
IntegrateTable[x_^3/(a_ - x_^2)^(3/2), x_] /; FreeQ[a, x] := 2 Sqrt[a - x^2] + x^2/Sqrt[a - x^2];

(* Formula 213 *)
IntegrateTable[1/(x_^3 Sqrt[a_ - x_^2]), x_] /; FreeQ[a, x] :=
  -Sqrt[a - x^2]/(2 a x^2) - 1/(2 a Sqrt[a]) Log[(Sqrt[a] + Sqrt[a - x^2])/x];

(* Formula 214 *)
IntegrateTable[1/(x_ (a_ - x_^2)^(3/2)), x_] /; FreeQ[a, x] :=
  1/(a Sqrt[a - x^2]) - 1/(a Sqrt[a]) Log[(Sqrt[a] + Sqrt[a - x^2])/x];

(* Formula 215 *)
IntegrateTable[1/(x_^2 (a_ - x_^2)^(3/2)), x_] /; FreeQ[a, x] :=
  1/a^2 (-Sqrt[a - x^2]/x + x/Sqrt[a - x^2]);

(* Formula 216.  (The CRC reference numerator/denominator and log argument are
   mis-transcribed; re-derived here from the antiderivative.) *)
IntegrateTable[1/(x_^3 (a_ - x_^2)^(3/2)), x_] /; FreeQ[a, x] :=
  (3 x^2 - a)/(2 a^2 x^2 Sqrt[a - x^2]) - 3/(4 a^2 Sqrt[a]) Log[(Sqrt[a] + Sqrt[a - x^2])/(Sqrt[a] - Sqrt[a - x^2])];

(* Formula 217 & 220: General Reduction Rules *)
IntegrateTable[x_^m_/Sqrt[a_ - x_^2], x_] /; FreeQ[{a, m}, x] && m =!= 0 && IntegerQ[m] && m > 1 :=
  -x^(m - 1)/m Sqrt[a - x^2] + ((m - 1) a)/m IntegrateTable[x^(m - 2)/Sqrt[a - x^2], x];
IntegrateTable[1/(x_^m_ Sqrt[a_ - x_^2]), x_] /; FreeQ[{a, m}, x] && m =!= 1 && IntegerQ[m] && m > 2 :=
  -Sqrt[a - x^2]/((m - 1) a x^(m - 1)) + (m - 2)/((m - 1) a) IntegrateTable[1/(x^(m - 2) Sqrt[a - x^2]), x];

(* Formula 223.  Both squared constants a_^2, b_^2 become linear a_, b_ (each
   positive), recovered via Sqrt on the RHS; a^2 > b^2 <-> a > b. *)
IntegrateTable[1/((b_ - x_^2) Sqrt[a_ - x_^2]), x_] /; FreeQ[{a, b}, x] && a > b :=
  1/(2 Sqrt[b] Sqrt[a - b]) Log[(Sqrt[b] Sqrt[a - x^2] + x Sqrt[a - b])/(Sqrt[b] Sqrt[a - x^2] - x Sqrt[a - b])];
IntegrateTable[1/((b_ - x_^2) Sqrt[a_ - x_^2]), x_] /; FreeQ[{a, b}, x] && b > a :=
  1/(Sqrt[b] Sqrt[b - a]) ArcTan[(x Sqrt[b - a])/(Sqrt[b] Sqrt[a - x^2])];

(* Formula 224 *)
IntegrateTable[1/((b_ + x_^2) Sqrt[a_ - x_^2]), x_] /; FreeQ[{a, b}, x] && Not[TrueQ[b < 0]] :=
  1/(Sqrt[b] Sqrt[a + b]) ArcTan[(x Sqrt[a + b])/(Sqrt[b] Sqrt[a - x^2])];

(* Formula 225 *)
IntegrateTable[Sqrt[a_ - x_^2]/(b_ + x_^2), x_] /; FreeQ[{a, b}, x] && Not[TrueQ[b < 0]] :=
  Sqrt[a + b]/Sqrt[b] ArcSin[(x Sqrt[a + b])/(Sqrt[a] Sqrt[x^2 + b])] - ArcSin[x/Sqrt[a]];

(* Formula 226 *)
IntegrateTable[1/Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c}, x] && c > 0 := 
  1/Sqrt[c] Log[2 Sqrt[c (a + b x + c x^2)] + 2 c x + b];
IntegrateTable[1/Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c}, x] && c < 0 := 
  -1/Sqrt[-c] ArcSin[(2 c x + b)/Sqrt[b^2 - 4 a c]];

(* Formula 227 & 228 *)
IntegrateTable[1/(a_ + b_. x_ + c_. x_^2)^(3/2), x_] /; FreeQ[{a, b, c}, x] :=
  With[{X = a + b x + c x^2, q = 4 a c - b^2}, (2 (2 c x + b))/(q Sqrt[X])];
IntegrateTable[1/(a_ + b_. x_ + c_. x_^2)^(5/2), x_] /; FreeQ[{a, b, c}, x] :=
  With[{q = 4 a c - b^2},
    With[{X = a + b x + c x^2, k = (4 c)/q},
      (2 (2 c x + b))/(3 q Sqrt[X]) (1/X + 2 k)
    ]
  ];

(* Formula 229: Generic reduction for n (mapped from texts n) *)
IntegrateTable[1/((a_ + b_. x_ + c_. x_^2)^n_ Sqrt[a_ + b_. x_ + c_. x_^2]), x_] /; FreeQ[{a, b, c, n}, x] && n =!= 0 && IntegerQ[n] && n > 0 :=
  With[{q = 4 a c - b^2},
    With[{X = a + b x + c x^2, k = (4 c)/q},
      (2 (2 c x + b) Sqrt[X])/((2 n + 1) q X^(n + 1)) + (2 k n)/(2 n + 1) IntegrateTable[1/(X^(n - 1) Sqrt[X]), x]
    ]
  ];

(* Formula 230 - 232.
 * `k = (4 c)/q` references q, so it must be in an inner With over q --
 * Mathilda/Mathematica's With is simultaneous-binding, so a single
 * With[{q = 4 a c - b^2, k = (4 c)/q}, ...] would leave k bound to
 * (4 c)/q with q as a free global symbol. *)
IntegrateTable[Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c}, x] :=
  With[{q = 4 a c - b^2},
    With[{X = a + b x + c x^2, k = (4 c)/q},
      ((2 c x + b) Sqrt[X])/(4 c) + 1/(2 k) IntegrateTable[1/Sqrt[X], x]
    ]
  ];
IntegrateTable[(a_ + b_. x_ + c_. x_^2)^(3/2), x_] /; FreeQ[{a, b, c}, x] :=
  With[{q = 4 a c - b^2},
    With[{X = a + b x + c x^2, k = (4 c)/q},
      ((2 c x + b) Sqrt[X])/(8 c) (X + 3/(2 k)) + 3/(8 k^2) IntegrateTable[1/Sqrt[X], x]
    ]
  ];
IntegrateTable[(a_ + b_. x_ + c_. x_^2)^(5/2), x_] /; FreeQ[{a, b, c}, x] :=
  With[{q = 4 a c - b^2},
    With[{X = a + b x + c x^2, k = (4 c)/q},
      ((2 c x + b) Sqrt[X])/(12 c) (X^2 + (5 X)/(4 k) + 15/(8 k^2)) + 5/(16 k^3) IntegrateTable[1/Sqrt[X], x]
    ]
  ];

(* Formula 233: General reduction *)
IntegrateTable[(a_ + b_. x_ + c_. x_^2)^n_ Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c, n}, x] && IntegerQ[n] && n > 0 :=
  With[{q = 4 a c - b^2},
    With[{X = a + b x + c x^2, k = (4 c)/q},
      ((2 c x + b) X^n Sqrt[X])/(4 (n + 1) c) + (2 n + 1)/(2 (n + 1) k) IntegrateTable[X^(n - 1) Sqrt[X], x]
    ]
  ];

(* Formula 234 *)
IntegrateTable[x_/Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2}, 
    Sqrt[X]/c - b/(2 c) IntegrateTable[1/Sqrt[X], x]
  ];

(* Formula 235 & 236 *)
IntegrateTable[x_/(a_ + b_. x_ + c_. x_^2)^(3/2), x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2, q = 4 a c - b^2}, -(2 (b x + 2 a))/(q Sqrt[X])];
IntegrateTable[x_/((a_ + b_. x_ + c_. x_^2)^n_ Sqrt[a_ + b_. x_ + c_. x_^2]), x_] /; FreeQ[{a, b, c, n}, x] && n =!= 0 := 
  With[{X = a + b x + c x^2}, 
    -Sqrt[X]/((2 n - 1) c X^n) - b/(2 c) IntegrateTable[1/(X^n Sqrt[X]), x]
  ];

(* Formula 237 - 239 *)
IntegrateTable[x_^2/Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2}, 
    (x/(2 c) - (3 b)/(4 c^2)) Sqrt[X] + (3 b^2 - 4 a c)/(8 c^2) IntegrateTable[1/Sqrt[X], x]
  ];
IntegrateTable[x_^2/(a_ + b_. x_ + c_. x_^2)^(3/2), x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2, q = 4 a c - b^2}, 
    ((2 b^2 - 4 a c) x + 2 a b)/(c q Sqrt[X]) + 1/c IntegrateTable[1/Sqrt[X], x]
  ];
IntegrateTable[x_^2/((a_ + b_. x_ + c_. x_^2)^n_ Sqrt[a_ + b_. x_ + c_. x_^2]), x_] /; FreeQ[{a, b, c, n}, x] && n =!= 0 := 
  With[{X = a + b x + c x^2, q = 4 a c - b^2}, 
    ((2 b^2 - 4 a c) x + 2 a b)/((2 n - 1) c q X^(n - 1) Sqrt[X]) + (4 a c + (2 n - 3) b^2)/((2 n - 1) c q) IntegrateTable[1/(X^(n - 1) Sqrt[X]), x]
  ];

(* Formula 240 & 241 *)
IntegrateTable[x_^3/Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2}, 
    (x^2/(3 c) - (5 b x)/(12 c^2) + (5 b^2)/(8 c^3) - (2 a)/(3 c^2)) Sqrt[X] + ((3 a b)/(4 c^2) - (5 b^3)/(16 c^3)) IntegrateTable[1/Sqrt[X], x]
  ];
IntegrateTable[x_^n_/Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c, n}, x] && n =!= 0 && IntegerQ[n] && n > 0 :=
  With[{X = a + b x + c x^2},
    1/(n c) x^(n - 1) Sqrt[X] - ((2 n - 1) b)/(2 n c) IntegrateTable[x^(n - 1)/Sqrt[X], x] - ((n - 1) a)/(n c) IntegrateTable[x^(n - 2)/Sqrt[X], x]
  ];

(* Formula 242 - 245.
 * Integrate[x X^(3/2), x] = X^(5/2)/(5 c) - (b/(2 c)) Integrate[X^(3/2), x],
 * derived from x = (X' - b)/(2 c) so x X^(3/2) = X^(3/2) X'/(2 c) - b X^(3/2)/(2 c),
 * and the first term integrates to (2/5) X^(5/2)/(2 c).  Folding in Formula 231
 * for Integrate[X^(3/2), x] gives the closed form below.  Verified by direct
 * differentiation: D[F, x] - x X^(3/2) reduces to 0.
 *
 * The previous body
 *   (X Sqrt[X])/(3 c) - (b (2 c x + b) Sqrt[X])/(8 c^2) - b/(4 c k) Int[1/Sqrt[X], x]
 * is the closed form for Integrate[x Sqrt[X], x] (^(1/2), NOT ^(3/2)) -- the
 * pattern's (3/2) exponent and the body had drifted out of sync. *)
IntegrateTable[x_ (a_ + b_. x_ + c_. x_^2)^(3/2), x_] /; FreeQ[{a, b, c}, x] :=
  With[{q = 4 a c - b^2},
    With[{X = a + b x + c x^2, k = (4 c)/q},
      X^(5/2)/(5 c)
        - (b (2 c x + b) Sqrt[X])/(16 c^2) (X + 3/(2 k))
        - (3 b)/(16 c k^2) IntegrateTable[1/Sqrt[X], x]
    ]
  ];
IntegrateTable[x_ X_^n_ Sqrt[X_], x_] /; FreeQ[{a, b, c, n}, x] && IntegerQ[n] && n > 0 :=
  With[{X = a + b x + c x^2},
    (X^(n + 1) Sqrt[X])/((2 n + 3) c) - b/(2 c) IntegrateTable[X^n Sqrt[X], x]
  ];
IntegrateTable[x_^2 Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2}, 
    (x - (5 b)/(6 c)) (X Sqrt[X])/(4 c) + (5 b^2 - 4 a c)/(16 c^2) IntegrateTable[Sqrt[X], x]
  ];

(* Formula 246 *)
IntegrateTable[1/(x_ Sqrt[a_ + b_. x_ + c_. x_^2]), x_] /; FreeQ[{a, b, c}, x] && a > 0 := 
  With[{X = a + b x + c x^2}, -1/Sqrt[a] Log[(2 Sqrt[a X] + 2 a + b x)/x]];
IntegrateTable[1/(x_ Sqrt[a_ + b_. x_ + c_. x_^2]), x_] /; FreeQ[{a, b, c}, x] && a < 0 := 
  1/Sqrt[-a] ArcSin[(b x + 2 a)/(x Sqrt[b^2 - 4 a c])];

(* Formula 247 - 249 *)
IntegrateTable[1/(x_^2 Sqrt[a_ + b_. x_ + c_. x_^2]), x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2}, -Sqrt[X]/(a x) - b/(2 a) IntegrateTable[1/(x Sqrt[X]), x]];
IntegrateTable[Sqrt[a_ + b_. x_ + c_. x_^2]/x_, x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2}, Sqrt[X] + b/2 IntegrateTable[1/Sqrt[X], x] + a IntegrateTable[1/(x Sqrt[X]), x]];
IntegrateTable[Sqrt[a_ + b_. x_ + c_. x_^2]/x_^2, x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2}, -Sqrt[X]/x + b/2 IntegrateTable[1/(x Sqrt[X]), x] + c IntegrateTable[1/Sqrt[X], x]];

(* Formula 250 *)
IntegrateTable[Sqrt[2 a_ x_ - x_^2], x_] /; FreeQ[a, x] := 
  1/2 ((x - a) Sqrt[2 a x - x^2] + a^2 ArcSin[(x - a)/Abs[a]]);

(* Formula 251 *)
IntegrateTable[1/Sqrt[2 a_ x_ - x_^2], x_] /; FreeQ[a, x] := 
  ArcSin[(x - a)/Abs[a]];

(* Formula 252 *)
IntegrateTable[x_^n_ Sqrt[2 a_ x_ - x_^2], x_] /; FreeQ[{a, n}, x] && n =!= -2 && IntegerQ[n] && n > 0 :=
  -((x^(n - 1) (2 a x - x^2)^(3/2))/(n + 2)) + ((2 n + 1) a)/(n + 2) IntegrateTable[x^(n - 1) Sqrt[2 a x - x^2], x];

(* Formula 253 *)
IntegrateTable[Sqrt[2 a_ x_ - x_^2]/x_^n_, x_] /; FreeQ[{a, n}, x] && n =!= 3/2 && IntegerQ[n] && n > 0 :=
  (2 a x - x^2)^(3/2)/((3 - 2 n) a x^n) + (n - 3)/((2 n - 3) a) IntegrateTable[Sqrt[2 a x - x^2]/x^(n - 1), x];

(* Formula 254 *)
IntegrateTable[x_^n_/Sqrt[2 a_ x_ - x_^2], x_] /; FreeQ[{a, n}, x] && n =!= 0 && IntegerQ[n] && n > 0 :=
  -(x^(n - 1) Sqrt[2 a x - x^2])/n + (a (2 n - 1))/n IntegrateTable[x^(n - 1)/Sqrt[2 a x - x^2], x];

(* Formula 255 *)
IntegrateTable[1/(x_^n_ Sqrt[2 a_ x_ - x_^2]), x_] /; FreeQ[{a, n}, x] && n =!= 1/2 && IntegerQ[n] && n > 0 :=
  Sqrt[2 a x - x^2]/(a (1 - 2 n) x^n) + (n - 1)/((2 n - 1) a) IntegrateTable[1/(x^(n - 1) Sqrt[2 a x - x^2]), x];

(* Formula 256 *)
IntegrateTable[1/(2 a_ x_ - x_^2)^(3/2), x_] /; FreeQ[a, x] := 
  (x - a)/(a^2 Sqrt[2 a x - x^2]);

(* Formula 257 *)
IntegrateTable[x_/(2 a_ x_ - x_^2)^(3/2), x_] /; FreeQ[a, x] := 
  x/(a Sqrt[2 a x - x^2]);

(* Formula 258 *)
IntegrateTable[1/Sqrt[2 a_ x_ + x_^2], x_] /; FreeQ[a, x] := 
  Log[x + a + Sqrt[2 a x + x^2]];

(* Formula 259 *)
IntegrateTable[Sqrt[a_. x_^2 + c_], x_] /; FreeQ[{a, c}, x] && a < 0 := 
  x/2 Sqrt[a x^2 + c] + c/(2 Sqrt[-a]) ArcSin[x Sqrt[-a/c]];
IntegrateTable[Sqrt[a_. x_^2 + c_], x_] /; FreeQ[{a, c}, x] && a > 0 := 
  x/2 Sqrt[a x^2 + c] + c/(2 Sqrt[a]) Log[x Sqrt[a] + Sqrt[a x^2 + c]];

(* Formula 260: branch-correct.  The classical form ArcSin[x] - Sqrt[1 - x^2]
   is right on (-1, 1) but its derivative (1 + x)/Sqrt[1 - x^2] cannot be
   shown equal to Sqrt[(1 + x)/(1 - x)] without PowerExpand (the radicals
   live on different branches).  Keeping Sqrt[(1 + x)/(1 - x)] literally
   in the antiderivative makes the diff D[F, x] - integrand collapse to 0
   under ordinary Simplify. *)
IntegrateTable[Sqrt[(1 + x_)/(1 - x_)], x_] :=
  (x - 1) Sqrt[(1 + x)/(1 - x)] + ArcSin[x];

(* Formula 261 *)
IntegrateTable[1/(x_ Sqrt[a_. x_^n_ + c_]), x_] /; FreeQ[{a, c, n}, x] && c > 0 := 
  1/(n Sqrt[c]) Log[(Sqrt[a x^n + c] - Sqrt[c])/(Sqrt[a x^n + c] + Sqrt[c])];
IntegrateTable[1/(x_ Sqrt[a_. x_^n_ + c_]), x_] /; FreeQ[{a, c, n}, x] && c < 0 := 
  2/(n Sqrt[-c]) ArcSec[Sqrt[-(a x^n)/c]];

(* Formula 262 *)
IntegrateTable[1/Sqrt[a_. x_^2 + c_], x_] /; FreeQ[{a, c}, x] && a < 0 := 
  1/Sqrt[-a] ArcSin[x Sqrt[-a/c]];
IntegrateTable[1/Sqrt[a_. x_^2 + c_], x_] /; FreeQ[{a, c}, x] && a > 0 := 
  1/Sqrt[a] Log[x Sqrt[a] + Sqrt[a x^2 + c]];

(* Formula 263: Generalized m+1/2 exponent to general m_ *)
IntegrateTable[(a_. x_^2 + c_)^m_, x_] /; FreeQ[{a, c, m}, x] && m =!= -1/2 && IntegerQ[m] && m > 0 :=
  (x (a x^2 + c)^m)/(2 m + 1) + (2 m c)/(2 m + 1) IntegrateTable[(a x^2 + c)^(m - 1), x];

(* Formula 264 *)
IntegrateTable[x_ (a_. x_^2 + c_)^m_, x_] /; FreeQ[{a, c, m}, x] && m =!= -1 := 
  (a x^2 + c)^(m + 1)/(2 a (m + 1));

(* Formula 265 *)
IntegrateTable[(a_. x_^2 + c_)^m_/x_, x_] /; FreeQ[{a, c, m}, x] && m =!= 0 && IntegerQ[m] && m > 0 :=
  (a x^2 + c)^m/(2 m) + c IntegrateTable[(a x^2 + c)^(m - 1)/x, x];

(* Formula 266 *)
IntegrateTable[1/(a_. x_^2 + c_)^m_, x_] /; FreeQ[{a, c, m}, x] && m =!= 1 && IntegerQ[m] && m > 1 :=
  x/((2 m - 2) c (a x^2 + c)^(m - 1)) + (2 m - 3)/((2 m - 2) c) IntegrateTable[1/(a x^2 + c)^(m - 1), x];

(* Formula 267: Already covered structurally by 261 out of the box! *)

(* Formula 268 *)
IntegrateTable[(1 + x_^2)/((1 - x_^2) Sqrt[1 + x_^4]), x_] := 
  1/Sqrt[2] Log[(x Sqrt[2] + Sqrt[1 + x^4])/(1 - x^2)];

(* Formula 269 *)
IntegrateTable[(1 - x_^2)/((1 + x_^2) Sqrt[1 + x_^4]), x_] := 
  1/Sqrt[2] ArcTan[(x Sqrt[2])/Sqrt[1 + x^4]];

(* Formula 270.  Squared constant a_^2 -> linear a_ (Sqrt[a] on the RHS); the
   x^n - a^2 sibling (271) shares the x_^n_ + a_ form under TrueQ[a < 0]. *)
IntegrateTable[1/(x_ Sqrt[x_^n_ + a_]), x_] /; FreeQ[{a, n}, x] && Not[TrueQ[a < 0]] :=
  -2/(n Sqrt[a]) Log[(Sqrt[a] + Sqrt[x^n + a])/Sqrt[x^n]];

(* Formula 271 *)
IntegrateTable[1/(x_ Sqrt[x_^n_ + a_]), x_] /; FreeQ[{a, n}, x] && TrueQ[a < 0] :=
  -2/(n Sqrt[-a]) ArcSin[Sqrt[-a]/Sqrt[x^n]];

(* Formula 272.  Cubed constant a_^3 -> linear a_; (x/a_orig)^(3/2) = Sqrt[x^3/a]. *)
IntegrateTable[Sqrt[x_/(a_ - x_^3)], x_] /; FreeQ[a, x] :=
  2/3 ArcSin[Sqrt[x^3/a]];

(* Formula 273 *)
IntegrateTable[Sin[a_. x_], x_] /; FreeQ[a, x] := -Cos[a x]/a;

(* Formula 274 *)
IntegrateTable[Cos[a_. x_], x_] /; FreeQ[a, x] := Sin[a x]/a;

(* Formula 275 *)
IntegrateTable[Tan[a_. x_], x_] /; FreeQ[a, x] := -Log[Cos[a x]]/a;

(* Formula 276 *)
IntegrateTable[Cot[a_. x_], x_] /; FreeQ[a, x] := Log[Sin[a x]]/a;

(* Formula 277 *)
IntegrateTable[Sec[a_. x_], x_] /; FreeQ[a, x] := Log[Sec[a x] + Tan[a x]]/a;

(* Formula 278 *)
IntegrateTable[Csc[a_. x_], x_] /; FreeQ[a, x] := Log[Csc[a x] - Cot[a x]]/a;

(* Formula 279 *)
IntegrateTable[Sin[a_. x_]^2, x_] /; FreeQ[a, x] := x/2 - Sin[2 a x]/(4 a);

(* Formula 280 *)
IntegrateTable[Sin[a_. x_]^3, x_] /; FreeQ[a, x] := -1/(3 a) Cos[a x] (Sin[a x]^2 + 2);

(* Formula 281 *)
IntegrateTable[Sin[a_. x_]^4, x_] /; FreeQ[a, x] := (3 x)/8 - Sin[2 a x]/(4 a) + Sin[4 a x]/(32 a);

(* Formula 282 *)
IntegrateTable[Sin[a_. x_]^n_, x_] /; FreeQ[{a, n}, x] && n =!= 0 && IntegerQ[n] && n > 1 :=
  -(Sin[a x]^(n - 1) Cos[a x])/(n a) + (n - 1)/n IntegrateTable[Sin[a x]^(n - 2), x];

(* Formula 285 *)
IntegrateTable[Cos[a_. x_]^2, x_] /; FreeQ[a, x] := x/2 + Sin[2 a x]/(4 a);

(* Formula 286 *)
IntegrateTable[Cos[a_. x_]^3, x_] /; FreeQ[a, x] := 1/(3 a) Sin[a x] (Cos[a x]^2 + 2);

(* Formula 287 *)
IntegrateTable[Cos[a_. x_]^4, x_] /; FreeQ[a, x] := (3 x)/8 + Sin[2 a x]/(4 a) + Sin[4 a x]/(32 a);

(* Formula 288 *)
IntegrateTable[Cos[a_. x_]^n_, x_] /; FreeQ[{a, n}, x] && n =!= 0 && IntegerQ[n] && n > 1 :=
  (Cos[a x]^(n - 1) Sin[a x])/(n a) + (n - 1)/n IntegrateTable[Cos[a x]^(n - 2), x];

(* Formula 291 *)
IntegrateTable[1/Sin[a_. x_]^2, x_] /; FreeQ[a, x] := -Cot[a x]/a;

(* Formula 292 *)
IntegrateTable[1/Sin[a_. x_]^m_, x_] /; FreeQ[{a, m}, x] && m =!= 1 && IntegerQ[m] && m > 2 :=
  -Cos[a x]/(a (m - 1) Sin[a x]^(m - 1)) + (m - 2)/(m - 1) IntegrateTable[1/Sin[a x]^(m - 2), x];

(* Formula 295 *)
IntegrateTable[1/Cos[a_. x_]^2, x_] /; FreeQ[a, x] := Tan[a x]/a;

(* Formula 296 *)
IntegrateTable[1/Cos[a_. x_]^m_, x_] /; FreeQ[{a, m}, x] && m =!= 1 && IntegerQ[m] && m > 2 :=
  Sin[a x]/(a (m - 1) Cos[a x]^(m - 1)) + (m - 2)/(m - 1) IntegrateTable[1/Cos[a x]^(m - 2), x];

(* Sec Reduction (missing from CRC tables; the reciprocal-power rules
   above only fire on 1/Cos[a x]^m, not on the Sec head).  Mirrors the
   Sech reduction; terminates at Sec[a x] (odd n) or a constant (even n). *)
IntegrateTable[Sec[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 1 :=
  (Sec[a x]^(n - 2) Tan[a x])/(a (n - 1)) + (n - 2)/(n - 1) IntegrateTable[Sec[a x]^(n - 2), x];

(* Csc Reduction (missing from CRC tables; the reciprocal-power rules
   above only fire on 1/Sin[a x]^m, not on the Csc head).  The second
   term is + (n-2)/(n-1), unlike the Csch analogue: the circular
   identity is Csc^2 = 1 + Cot^2 whereas Csch^2 = Coth^2 - 1, which
   flips the recursion sign.  Terminates at Csc[a x] (odd n) or a
   constant (even n). *)
IntegrateTable[Csc[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 1 :=
  -(Csc[a x]^(n - 2) Cot[a x])/(a (n - 1)) + (n - 2)/(n - 1) IntegrateTable[Csc[a x]^(n - 2), x];

(* Formula 299 *)
IntegrateTable[Sin[m_. x_] Sin[n_. x_], x_] /; FreeQ[{m, n}, x] && m^2 =!= n^2 := 
  Sin[(m - n) x]/(2 (m - n)) - Sin[(m + n) x]/(2 (m + n));

(* Formula 300 *)
IntegrateTable[Cos[m_. x_] Cos[n_. x_], x_] /; FreeQ[{m, n}, x] && m^2 =!= n^2 := 
  Sin[(m - n) x]/(2 (m - n)) + Sin[(m + n) x]/(2 (m + n));

(* Formula 301 *)
IntegrateTable[Sin[a_. x_] Cos[a_. x_], x_] /; FreeQ[a, x] := Sin[a x]^2/(2 a);

(* Formula 302 *)
IntegrateTable[Sin[m_. x_] Cos[n_. x_], x_] /; FreeQ[{m, n}, x] && m^2 =!= n^2 := 
  -Cos[(m - n) x]/(2 (m - n)) - Cos[(m + n) x]/(2 (m + n));

(* Formula 303 *)
IntegrateTable[Sin[a_. x_]^2 Cos[a_. x_]^2, x_] /; FreeQ[a, x] := x/8 - Sin[4 a x]/(32 a);

(* Formula 304 *)
IntegrateTable[Sin[a_. x_] Cos[a_. x_]^m_, x_] /; FreeQ[{a, m}, x] && m =!= -1 := 
  -Cos[a x]^(m + 1)/((m + 1) a);

(* Formula 305 *)
IntegrateTable[Sin[a_. x_]^m_ Cos[a_. x_], x_] /; FreeQ[{a, m}, x] && m =!= -1 := 
  Sin[a x]^(m + 1)/((m + 1) a);

(* Formula 306 *)
IntegrateTable[Cos[a_. x_]^m_ Sin[a_. x_]^n_, x_] /; FreeQ[{a, m, n}, x] && m + n =!= 0 && IntegerQ[m] && m > 1 :=
  (Cos[a x]^(m - 1) Sin[a x]^(n + 1))/((m + n) a) + (m - 1)/(m + n) IntegrateTable[Cos[a x]^(m - 2) Sin[a x]^n, x];

(* Formula 307 *)
IntegrateTable[Cos[a_. x_]^m_/Sin[a_. x_]^n_, x_] /; FreeQ[{a, m, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 :=
  -Cos[a x]^(m + 1)/(a (n - 1) Sin[a x]^(n - 1)) - (m - n + 2)/(n - 1) IntegrateTable[Cos[a x]^m/Sin[a x]^(n - 2), x];

(* Formula 308 *)
IntegrateTable[Sin[a_. x_]^m_/Cos[a_. x_]^n_, x_] /; FreeQ[{a, m, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 :=
  Sin[a x]^(m - 1)/(a (n - 1) Cos[a x]^(n - 1)) - (m - n + 2)/(n - 1) IntegrateTable[Sin[a x]^m/Cos[a x]^(n - 2), x];

(* Mixed tangent/secant and cotangent/cosecant powers.  Mathilda canonicalises
   the Formula 307/308 quotients into these heads: Sin^m/Cos^n -> Tan^m Sec^(n-m)
   and Cos^m/Sin^n -> Cot^m Csc^(n-m), so the quotient rules above never match an
   evaluated integrand.  Reduce the Tan/Cot power by two (bottoming out on the
   single-Sec^n / Csc^n rules), then close the odd tail on the Tan Sec^n /
   Cot Csc^n bases (which include Sec Tan and Csc Cot at n = 1 — a form the CRC
   circular tables omit but the hyperbolic Sech Tanh / Csch Coth carry). *)
IntegrateTable[Tan[a_. x_] Sec[a_. x_]^n_., x_] /; FreeQ[{a, n}, x] && IntegerQ[n] :=
  Sec[a x]^n/(a n);
IntegrateTable[Cot[a_. x_] Csc[a_. x_]^n_., x_] /; FreeQ[{a, n}, x] && IntegerQ[n] :=
  -Csc[a x]^n/(a n);
IntegrateTable[Tan[a_. x_]^m_Integer Sec[a_. x_]^n_., x_] /; FreeQ[{a, m, n}, x] && IntegerQ[n] && m > 1 :=
  Tan[a x]^(m - 1) Sec[a x]^n/(a (m + n - 1)) - (m - 1)/(m + n - 1) IntegrateTable[Tan[a x]^(m - 2) Sec[a x]^n, x];
IntegrateTable[Cot[a_. x_]^m_Integer Csc[a_. x_]^n_., x_] /; FreeQ[{a, m, n}, x] && IntegerQ[n] && m > 1 :=
  -Cot[a x]^(m - 1) Csc[a x]^n/(a (m + n - 1)) - (m - 1)/(m + n - 1) IntegrateTable[Cot[a x]^(m - 2) Csc[a x]^n, x];

(* Formula 309 *)
IntegrateTable[Sin[a_. x_]/Cos[a_. x_]^2, x_] /; FreeQ[a, x] := Sec[a x]/a;

(* Formula 310 *)
IntegrateTable[Sin[a_. x_]^2/Cos[a_. x_], x_] /; FreeQ[a, x] := 
  -Sin[a x]/a + 1/a Log[Tan[Pi/4 + (a x)/2]];

(* Formula 311 *)
IntegrateTable[Cos[a_. x_]/Sin[a_. x_]^2, x_] /; FreeQ[a, x] := -Csc[a x]/a;

(* Formula 312 *)
IntegrateTable[1/(Sin[a_. x_] Cos[a_. x_]), x_] /; FreeQ[a, x] := Log[Tan[a x]]/a;

(* Formula 313 *)
IntegrateTable[1/(Sin[a_. x_] Cos[a_. x_]^2), x_] /; FreeQ[a, x] := 
  1/a (Sec[a x] + Log[Tan[(a x)/2]]);

(* Formula 314 *)
IntegrateTable[1/(Sin[a_. x_] Cos[a_. x_]^n_), x_] /; FreeQ[{a, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 :=
  1/(a (n - 1) Cos[a x]^(n - 1)) + IntegrateTable[1/(Sin[a x] Cos[a x]^(n - 2)), x];

(* Formula 315 *)
IntegrateTable[1/(Sin[a_. x_]^2 Cos[a_. x_]), x_] /; FreeQ[a, x] := 
  -Csc[a x]/a + 1/a Log[Tan[Pi/4 + (a x)/2]];

(* Formula 316 *)
IntegrateTable[1/(Sin[a_. x_]^2 Cos[a_. x_]^2), x_] /; FreeQ[a, x] := -2/a Cot[2 a x];

(* Formula 317 *)
IntegrateTable[1/(Sin[a_. x_]^m_ Cos[a_. x_]^n_), x_] /; FreeQ[{a, m, n}, x] && m =!= 1 && IntegerQ[m] && m > 1 :=
  -(1/(a (m - 1) Sin[a x]^(m - 1) Cos[a x]^(n - 1))) + (m + n - 2)/(m - 1) IntegrateTable[1/(Sin[a x]^(m - 2) Cos[a x]^n), x];

(* Formulas 312-317 keyed on the Csc / Sec HEADS.  Mathilda canonicalises
   1/Sin -> Csc and 1/Cos -> Sec, so the reciprocal `1/(Sin^m Cos^n)` forms
   above never match an evaluated integrand (e.g. Csc[x]^4 Sec[x]^3); these
   head-basis reductions supply the coverage.  The recursive term is
   + (m+n-2)/(m-1) (from 1 = Sin^2 + Cos^2), the opposite sign to the
   hyperbolic Csch/Sech analogue.  Terminate at Csc Sec or single Sec^n/Csc^n. *)
IntegrateTable[Csc[a_. x_] Sec[a_. x_], x_] /; FreeQ[a, x] := Log[Tan[a x]]/a;
IntegrateTable[Csc[a_. x_] Sec[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 1 :=
  Sec[a x]^(n - 1)/(a (n - 1)) + IntegrateTable[Csc[a x] Sec[a x]^(n - 2), x];
IntegrateTable[Csc[a_. x_]^m_Integer Sec[a_. x_]^n_Integer, x_] /; FreeQ[{a, m, n}, x] && m > 1 :=
  -(Csc[a x]^(m - 1) Sec[a x]^(n - 1))/(a (m - 1)) + (m + n - 2)/(m - 1) IntegrateTable[Csc[a x]^(m - 2) Sec[a x]^n, x];

(* Formula 318 *)
IntegrateTable[Sin[a_ + b_. x_], x_] /; FreeQ[{a, b}, x] := -Cos[a + b x]/b;

(* Formula 319 *)
IntegrateTable[Cos[a_ + b_. x_], x_] /; FreeQ[{a, b}, x] := Sin[a + b x]/b;

(* Formula 320 *)
IntegrateTable[1/(1 + Sin[a_. x_]), x_] /; FreeQ[a, x] := -1/a Tan[Pi/4 - (a x)/2];
IntegrateTable[1/(1 - Sin[a_. x_]), x_] /; FreeQ[a, x] := 1/a Tan[Pi/4 + (a x)/2];

(* Formula 321 *)
IntegrateTable[1/(1 + Cos[a_. x_]), x_] /; FreeQ[a, x] := 
  1/a Tan[(a x)/2];

(* Formula 322 *)
IntegrateTable[1/(1 - Cos[a_. x_]), x_] /; FreeQ[a, x] := 
  -1/a Cot[(a x)/2];

(* Formula 323 *)
IntegrateTable[1/(a_ + b_. Sin[x_]), x_] /; FreeQ[{a, b}, x] && a^2 > b^2 := 
  2/Sqrt[a^2 - b^2] ArcTan[(a Tan[x/2] + b)/Sqrt[a^2 - b^2]];
IntegrateTable[1/(a_ + b_. Sin[x_]), x_] /; FreeQ[{a, b}, x] && b^2 > a^2 := 
  1/Sqrt[b^2 - a^2] Log[(a Tan[x/2] + b - Sqrt[b^2 - a^2])/(a Tan[x/2] + b + Sqrt[b^2 - a^2])];

(* Formula 324 *)
IntegrateTable[1/(a_ + b_. Cos[x_]), x_] /; FreeQ[{a, b}, x] && a^2 > b^2 := 
  2/Sqrt[a^2 - b^2] ArcTan[(Sqrt[a^2 - b^2] Tan[x/2])/(a + b)];
IntegrateTable[1/(a_ + b_. Cos[x_]), x_] /; FreeQ[{a, b}, x] && b^2 > a^2 := 
  1/Sqrt[b^2 - a^2] Log[(Sqrt[b^2 - a^2] Tan[x/2] + a + b)/(Sqrt[b^2 - a^2] Tan[x/2] - a - b)];

(* Formula 327 *)
IntegrateTable[1/(a_ Cos[x_]^2 + b_ Sin[x_]^2), x_] /; FreeQ[{a, b}, x] :=
  1/(Sqrt[a] Sqrt[b]) ArcTan[(Sqrt[b] Tan[x])/Sqrt[a]];

(* Formula 329 *)
IntegrateTable[(Sin[c_. x_] Cos[c_. x_])/(a_ Cos[c_. x_]^2 + b_ Sin[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] && a =!= b := 
  1/(2 c (b - a)) Log[a Cos[c x]^2 + b Sin[c x]^2];

(* Formula 330 *)
IntegrateTable[1/(a_ + b_. Tan[c_. x_]), x_] /; FreeQ[{a, b, c}, x] := 
  1/(c (a^2 + b^2)) (a c x + b Log[a Cos[c x] + b Sin[c x]]);

(* Formula 331 *)
IntegrateTable[1/(b_ + a_. Cot[c_. x_]), x_] /; FreeQ[{a, b, c}, x] := 
  1/(c (a^2 + b^2)) (b c x - a Log[a Cos[c x] + b Sin[c x]]);

(* Formula 333 *)
IntegrateTable[Sin[a_. x_]/(1 + Sin[a_. x_]), x_] /; FreeQ[a, x] := 
  x + 1/a Tan[Pi/4 - (a x)/2];
IntegrateTable[Sin[a_. x_]/(1 - Sin[a_. x_]), x_] /; FreeQ[a, x] := 
  -x + 1/a Tan[Pi/4 + (a x)/2];

(* Formula 334 *)
IntegrateTable[1/(Sin[a_. x_] (1 + Sin[a_. x_])), x_] /; FreeQ[a, x] := 
  1/a Tan[Pi/4 - (a x)/2] + 1/a Log[Tan[(a x)/2]];
IntegrateTable[1/(Sin[a_. x_] (1 - Sin[a_. x_])), x_] /; FreeQ[a, x] := 
  1/a Tan[Pi/4 + (a x)/2] + 1/a Log[Tan[(a x)/2]];

(* Formula 335 & 336 *)
IntegrateTable[1/(1 + Sin[a_. x_])^2, x_] /; FreeQ[a, x] := 
  -1/(2 a) Tan[Pi/4 - (a x)/2] - 1/(6 a) Tan[Pi/4 - (a x)/2]^3;
IntegrateTable[1/(1 - Sin[a_. x_])^2, x_] /; FreeQ[a, x] := 
  1/(2 a) Cot[Pi/4 - (a x)/2] + 1/(6 a) Cot[Pi/4 - (a x)/2]^3;

(* Formula 337 & 338 *)
IntegrateTable[Sin[a_. x_]/(1 + Sin[a_. x_])^2, x_] /; FreeQ[a, x] := 
  -1/(2 a) Tan[Pi/4 - (a x)/2] + 1/(6 a) Tan[Pi/4 - (a x)/2]^3;
IntegrateTable[Sin[a_. x_]/(1 - Sin[a_. x_])^2, x_] /; FreeQ[a, x] := 
  -1/(2 a) Cot[Pi/4 - (a x)/2] + 1/(6 a) Cot[Pi/4 - (a x)/2]^3;

(* Formula 339 *)
IntegrateTable[Sin[x_]/(a_ + b_. Sin[x_]), x_] /; FreeQ[{a, b}, x] := 
  x/b - a/b IntegrateTable[1/(a + b Sin[x]), x];

(* Formula 340 *)
IntegrateTable[1/(Sin[x_] (a_ + b_. Sin[x_])), x_] /; FreeQ[{a, b}, x] := 
  1/a Log[Tan[x/2]] - b/a IntegrateTable[1/(a + b Sin[x]), x];

(* Formula 341 *)
IntegrateTable[1/(a_ + b_. Sin[x_])^2, x_] /; FreeQ[{a, b}, x] && a^2 =!= b^2 := 
  (b Cos[x])/((a^2 - b^2) (a + b Sin[x])) + a/(a^2 - b^2) IntegrateTable[1/(a + b Sin[x]), x];

(* Formula 342 & 343 *)
IntegrateTable[1/(a_ + b_ Sin[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] && Not[TrueQ[b < 0]] :=
  1/(Sqrt[a] c Sqrt[a + b]) ArcTan[(Sqrt[a + b] Tan[c x])/Sqrt[a]];
IntegrateTable[1/(a_ + b_ Sin[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] && TrueQ[b < 0] && a + b > 0 :=
  1/(Sqrt[a] c Sqrt[a + b]) ArcTan[(Sqrt[a + b] Tan[c x])/Sqrt[a]];
IntegrateTable[1/(a_ + b_ Sin[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] && TrueQ[b < 0] && a + b < 0 :=
  1/(2 Sqrt[a] c Sqrt[-a - b]) Log[(Sqrt[-a - b] Tan[c x] + Sqrt[a])/(Sqrt[-a - b] Tan[c x] - Sqrt[a])];

(* Formula 344 & 345 *)
IntegrateTable[Cos[a_. x_]/(1 + Cos[a_. x_]), x_] /; FreeQ[a, x] := 
  x - 1/a Tan[(a x)/2];
IntegrateTable[Cos[a_. x_]/(1 - Cos[a_. x_]), x_] /; FreeQ[a, x] := 
  -x - 1/a Cot[(a x)/2];

(* Formula 346 & 347 *)
IntegrateTable[1/(Cos[a_. x_] (1 + Cos[a_. x_])), x_] /; FreeQ[a, x] := 
  1/a Log[Tan[Pi/4 + (a x)/2]] - 1/a Tan[(a x)/2];
IntegrateTable[1/(Cos[a_. x_] (1 - Cos[a_. x_])), x_] /; FreeQ[a, x] := 
  1/a Log[Tan[Pi/4 + (a x)/2]] - 1/a Cot[(a x)/2];

(* Formula 348 & 349 *)
IntegrateTable[1/(1 + Cos[a_. x_])^2, x_] /; FreeQ[a, x] := 
  1/(2 a) Tan[(a x)/2] + 1/(6 a) Tan[(a x)/2]^3;
IntegrateTable[1/(1 - Cos[a_. x_])^2, x_] /; FreeQ[a, x] := 
  -1/(2 a) Cot[(a x)/2] - 1/(6 a) Cot[(a x)/2]^3;

(* Formula 350 & 351 *)
IntegrateTable[Cos[a_. x_]/(1 + Cos[a_. x_])^2, x_] /; FreeQ[a, x] := 
  1/(2 a) Tan[(a x)/2] - 1/(6 a) Tan[(a x)/2]^3;
IntegrateTable[Cos[a_. x_]/(1 - Cos[a_. x_])^2, x_] /; FreeQ[a, x] := 
  1/(2 a) Cot[(a x)/2] - 1/(6 a) Cot[(a x)/2]^3;

(* Formula 352 *)
IntegrateTable[Cos[x_]/(a_ + b_. Cos[x_]), x_] /; FreeQ[{a, b}, x] := 
  x/b - a/b IntegrateTable[1/(a + b Cos[x]), x];

(* Formula 353 *)
IntegrateTable[1/(Cos[x_] (a_ + b_. Cos[x_])), x_] /; FreeQ[{a, b}, x] := 
  1/a Log[Tan[x/2 + Pi/4]] - b/a IntegrateTable[1/(a + b Cos[x]), x];

(* Formula 354 *)
IntegrateTable[1/(a_ + b_. Cos[x_])^2, x_] /; FreeQ[{a, b}, x] && a^2 =!= b^2 := 
  (b Sin[x])/((b^2 - a^2) (a + b Cos[x])) - a/(b^2 - a^2) IntegrateTable[1/(a + b Cos[x]), x];

(* Formula 355 *)
IntegrateTable[Cos[x_]/(a_ + b_. Cos[x_])^2, x_] /; FreeQ[{a, b}, x] && a^2 =!= b^2 := 
  (a Sin[x])/((a^2 - b^2) (a + b Cos[x])) - b/(a^2 - b^2) IntegrateTable[1/(a + b Cos[x]), x];

(* Formula 356 (1/(a_^2 + b_^2 - 2 a_ b_ Cos[c x])) omitted: the squared/product
   pattern constants never matched a numeric argument, and once expanded the
   integrand is an ordinary 1/(p + q Cos[c x]) form already handled above. *)

(* Formula 357 & 358 *)
IntegrateTable[1/(a_ + b_ Cos[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] && Not[TrueQ[b < 0]] :=
  1/(Sqrt[a] c Sqrt[a + b]) ArcTan[(Sqrt[a] Tan[c x])/Sqrt[a + b]];
IntegrateTable[1/(a_ + b_ Cos[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] && TrueQ[b < 0] && a + b > 0 :=
  1/(Sqrt[a] c Sqrt[a + b]) ArcTan[(Sqrt[a] Tan[c x])/Sqrt[a + b]];
IntegrateTable[1/(a_ + b_ Cos[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] && TrueQ[b < 0] && a + b < 0 :=
  1/(2 Sqrt[a] c Sqrt[-a - b]) Log[(Sqrt[a] Tan[c x] - Sqrt[-a - b])/(Sqrt[a] Tan[c x] + Sqrt[-a - b])];

(* Formula 359 *)
IntegrateTable[Sin[a_. x_]/(1 + Cos[a_. x_]), x_] /; FreeQ[a, x] := -1/a Log[1 + Cos[a x]];
IntegrateTable[Sin[a_. x_]/(1 - Cos[a_. x_]), x_] /; FreeQ[a, x] := 1/a Log[1 - Cos[a x]];

(* Formula 360 *)
IntegrateTable[Cos[a_. x_]/(1 + Sin[a_. x_]), x_] /; FreeQ[a, x] := 1/a Log[1 + Sin[a x]];
IntegrateTable[Cos[a_. x_]/(1 - Sin[a_. x_]), x_] /; FreeQ[a, x] := -1/a Log[1 - Sin[a x]];

(* Formula 361 *)
IntegrateTable[1/(Sin[a_. x_] (1 + Cos[a_. x_])), x_] /; FreeQ[a, x] := 
  1/(2 a (1 + Cos[a x])) + 1/(2 a) Log[Tan[(a x)/2]];
IntegrateTable[1/(Sin[a_. x_] (1 - Cos[a_. x_])), x_] /; FreeQ[a, x] := 
  -1/(2 a (1 - Cos[a x])) + 1/(2 a) Log[Tan[(a x)/2]];

(* Formula 362 *)
IntegrateTable[1/(Cos[a_. x_] (1 + Sin[a_. x_])), x_] /; FreeQ[a, x] := 
  -1/(2 a (1 + Sin[a x])) + 1/(2 a) Log[Tan[(a x)/2 + Pi/4]];
IntegrateTable[1/(Cos[a_. x_] (1 - Sin[a_. x_])), x_] /; FreeQ[a, x] := 
  1/(2 a (1 - Sin[a x])) + 1/(2 a) Log[Tan[(a x)/2 + Pi/4]];

(* Formula 363 *)
IntegrateTable[Sin[a_. x_]/(Cos[a_. x_] (1 + Cos[a_. x_])), x_] /; FreeQ[a, x] := 1/a Log[Sec[a x] + 1];
IntegrateTable[Sin[a_. x_]/(Cos[a_. x_] (1 - Cos[a_. x_])), x_] /; FreeQ[a, x] := 1/a Log[Sec[a x] - 1];

(* Formula 364 *)
IntegrateTable[Cos[a_. x_]/(Sin[a_. x_] (1 + Sin[a_. x_])), x_] /; FreeQ[a, x] := -1/a Log[Csc[a x] + 1];
IntegrateTable[Cos[a_. x_]/(Sin[a_. x_] (1 - Sin[a_. x_])), x_] /; FreeQ[a, x] := -1/a Log[Csc[a x] - 1];

(* Formula 365 *)
IntegrateTable[Sin[a_. x_]/(Cos[a_. x_] (1 + Sin[a_. x_])), x_] /; FreeQ[a, x] := 
  1/(2 a (1 + Sin[a x])) + 1/(2 a) Log[Tan[(a x)/2 + Pi/4]];
IntegrateTable[Sin[a_. x_]/(Cos[a_. x_] (1 - Sin[a_. x_])), x_] /; FreeQ[a, x] := 
  1/(2 a (1 - Sin[a x])) - 1/(2 a) Log[Tan[(a x)/2 + Pi/4]];

(* Formula 366 *)
IntegrateTable[Cos[a_. x_]/(Sin[a_. x_] (1 + Cos[a_. x_])), x_] /; FreeQ[a, x] := 
  -1/(2 a (1 + Cos[a x])) + 1/(2 a) Log[Tan[(a x)/2]];
IntegrateTable[Cos[a_. x_]/(Sin[a_. x_] (1 - Cos[a_. x_])), x_] /; FreeQ[a, x] := 
  -1/(2 a (1 - Cos[a x])) - 1/(2 a) Log[Tan[(a x)/2]];

(* Formula 367 *)
IntegrateTable[1/(Sin[a_. x_] + Cos[a_. x_]), x_] /; FreeQ[a, x] := 
  1/(a Sqrt[2]) Log[Tan[(a x)/2 + Pi/8]];
IntegrateTable[1/(Sin[a_. x_] - Cos[a_. x_]), x_] /; FreeQ[a, x] := 
  1/(a Sqrt[2]) Log[Tan[(a x)/2 - Pi/8]];

(* Formula 368 *)
IntegrateTable[1/(Sin[a_. x_] + Cos[a_. x_])^2, x_] /; FreeQ[a, x] := 
  1/(2 a) Tan[a x - Pi/4];
IntegrateTable[1/(Sin[a_. x_] - Cos[a_. x_])^2, x_] /; FreeQ[a, x] := 
  1/(2 a) Tan[a x + Pi/4];

(* Formula 369 *)
IntegrateTable[1/(1 + Cos[a_. x_] + Sin[a_. x_]), x_] /; FreeQ[a, x] := 
  1/a Log[1 + Tan[(a x)/2]];
IntegrateTable[1/(1 + Cos[a_. x_] - Sin[a_. x_]), x_] /; FreeQ[a, x] := 
  -1/a Log[1 - Tan[(a x)/2]];

(* Formula 370 *)
IntegrateTable[1/(a_ Cos[c_. x_]^2 + b_ Sin[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] && TrueQ[b < 0] :=
  1/(2 Sqrt[a] Sqrt[-b] c) Log[(Sqrt[-b] Tan[c x] + Sqrt[a])/(Sqrt[-b] Tan[c x] - Sqrt[a])];

(* Formula 371 *)
IntegrateTable[x_ Sin[a_. x_], x_] /; FreeQ[a, x] := 
  1/a^2 Sin[a x] - x/a Cos[a x];

(* Formula 372 *)
IntegrateTable[x_^2 Sin[a_. x_], x_] /; FreeQ[a, x] := 
  (2 x)/a^2 Sin[a x] + (2 - a^2 x^2)/a^3 Cos[a x];

(* Formula 373 *)
IntegrateTable[x_^3 Sin[a_. x_], x_] /; FreeQ[a, x] := 
  (3 a^2 x^2 - 6)/a^4 Sin[a x] + (6 x - a^2 x^3)/a^3 Cos[a x];

(* Formula 374: General Reduction *)
IntegrateTable[x_^m_ Sin[a_. x_], x_] /; FreeQ[{a, m}, x] && m > 0 := 
  -(x^m Cos[a x]/a) + m/a IntegrateTable[x^(m - 1) Cos[a x], x];

(* Formula 375 *)
IntegrateTable[x_ Cos[a_. x_], x_] /; FreeQ[a, x] := 
  1/a^2 Cos[a x] + x/a Sin[a x];

(* Formula 376 *)
IntegrateTable[x_^2 Cos[a_. x_], x_] /; FreeQ[a, x] := 
  (2 x)/a^2 Cos[a x] + (a^2 x^2 - 2)/a^3 Sin[a x];

(* Formula 377 *)
IntegrateTable[x_^3 Cos[a_. x_], x_] /; FreeQ[a, x] := 
  (3 a^2 x^2 - 6)/a^4 Cos[a x] + (a^2 x^3 - 6 x)/a^3 Sin[a x];

(* Formula 378: General Reduction *)
IntegrateTable[x_^m_ Cos[a_. x_], x_] /; FreeQ[{a, m}, x] && m > 0 := 
  x^m Sin[a x]/a - m/a IntegrateTable[x^(m - 1) Sin[a x], x];

(* Formula 381 *)
IntegrateTable[x_ Sin[a_. x_]^2, x_] /; FreeQ[a, x] := 
  x^2/4 - x/(4 a) Sin[2 a x] - 1/(8 a^2) Cos[2 a x];

(* Formula 382 *)
IntegrateTable[x_^2 Sin[a_. x_]^2, x_] /; FreeQ[a, x] := 
  x^3/6 - (x^2/(4 a) - 1/(8 a^3)) Sin[2 a x] - x/(4 a^2) Cos[2 a x];

(* Formula 383 *)
IntegrateTable[x_ Sin[a_. x_]^3, x_] /; FreeQ[a, x] := 
  x/(12 a) Cos[3 a x] - 1/(36 a^2) Sin[3 a x] - (3 x)/(4 a) Cos[a x] + 3/(4 a^2) Sin[a x];

(* Formula 384 *)
IntegrateTable[x_ Cos[a_. x_]^2, x_] /; FreeQ[a, x] := 
  x^2/4 + x/(4 a) Sin[2 a x] + 1/(8 a^2) Cos[2 a x];

(* Formula 385 *)
IntegrateTable[x_^2 Cos[a_. x_]^2, x_] /; FreeQ[a, x] := 
  x^3/6 + (x^2/(4 a) - 1/(8 a^3)) Sin[2 a x] + x/(4 a^2) Cos[2 a x];

(* Formula 386 *)
IntegrateTable[x_ Cos[a_. x_]^3, x_] /; FreeQ[a, x] := 
  x/(12 a) Sin[3 a x] + 1/(36 a^2) Cos[3 a x] + (3 x)/(4 a) Sin[a x] + 3/(4 a^2) Cos[a x];

(* Formula 387 base case (m = 1): the sine / cosine integrals.  These also
   terminate the Formula 387/388 negative-power recurrences below, which reduce
   Sin[a x]/x^m and Cos[a x]/x^m down to Sin[a x]/x and Cos[a x]/x.
   D[SinIntegral[a x], x] = Sin[a x]/x, D[CosIntegral[a x], x] = Cos[a x]/x. *)
IntegrateTable[Sin[a_. x_]/x_, x_] /; FreeQ[a, x] := SinIntegral[a x];
IntegrateTable[Cos[a_. x_]/x_, x_] /; FreeQ[a, x] := CosIntegral[a x];

(* Formula 387: Reduction for Negative Powers *)
IntegrateTable[Sin[a_. x_]/x_^m_, x_] /; FreeQ[{a, m}, x] && m =!= 1 && IntegerQ[m] && m > 1 :=
  -Sin[a x]/((m - 1) x^(m - 1)) + a/(m - 1) IntegrateTable[Cos[a x]/x^(m - 1), x];

(* Formula 388: Reduction for Negative Powers *)
IntegrateTable[Cos[a_. x_]/x_^m_, x_] /; FreeQ[{a, m}, x] && m =!= 1 && IntegerQ[m] && m > 1 :=
  -Cos[a x]/((m - 1) x^(m - 1)) - a/(m - 1) IntegrateTable[Sin[a x]/x^(m - 1), x];

(* Formula 389 *)
IntegrateTable[x_/(1 + Sin[a_. x_]), x_] /; FreeQ[a, x] := 
  -(x Cos[a x])/(a (1 + Sin[a x])) + 1/a^2 Log[1 + Sin[a x]];
IntegrateTable[x_/(1 - Sin[a_. x_]), x_] /; FreeQ[a, x] := 
  (x Cos[a x])/(a (1 - Sin[a x])) + 1/a^2 Log[1 - Sin[a x]];

(* Formula 390 *)
IntegrateTable[x_/(1 + Cos[a_. x_]), x_] /; FreeQ[a, x] := 
  x/a Tan[(a x)/2] + 2/a^2 Log[Cos[(a x)/2]];

(* Formula 391 *)
IntegrateTable[x_/(1 - Cos[a_. x_]), x_] /; FreeQ[a, x] := 
  -x/a Cot[(a x)/2] + 2/a^2 Log[Sin[(a x)/2]];

(* Formula 392 *)
IntegrateTable[(x_ + Sin[x_])/(1 + Cos[x_]), x_] := x Tan[x/2];

(* Formula 393 *)
IntegrateTable[(x_ - Sin[x_])/(1 - Cos[x_]), x_] := -x Cot[x/2];

(* Formula 394: branch-correct.  1 - Cos[a x] = 2 Sin[a x/2]^2, so the
   classical form (-2 Sqrt[2]/a) Cos[a x/2] differentiates back to
   Sqrt[2] |Sin[a x/2]|, matching the integrand only when Sin[a x/2] >= 0.
   The form below keeps the integrand's Sqrt[1 - Cos[a x]] literally,
   giving a derivative that equals the integrand on every branch. *)
IntegrateTable[Sqrt[1 - Cos[a_. x_]], x_] /; FreeQ[a, x] := -(2/a) Cot[(a x)/2] Sqrt[1 - Cos[a x]];

(* Formula 395: branch-correct counterpart.  See Formula 394. *)
IntegrateTable[Sqrt[1 + Cos[a_. x_]], x_] /; FreeQ[a, x] := (2/a) Tan[(a x)/2] Sqrt[1 + Cos[a x]];

(* Formula 396: branch-correct.  The classical primary-branch form
   2 (Sin[x/2] - Cos[x/2]) differentiates to Cos[x/2] + Sin[x/2], but
   Sqrt[1 + Sin[x]] = |Sin[x/2] + Cos[x/2]|, so it is correct only where
   that sum is >= 0.  Keeping the radical literally, -2 Cos[x]/Sqrt[1+Sin[x]]
   differentiates to (Sin^2 + 2 Sin + 1)/(1+Sin)^(3/2) = Sqrt[1+Sin[x]] on
   every branch. *)
IntegrateTable[Sqrt[1 + Sin[x_]], x_] := -((2 Cos[x])/Sqrt[1 + Sin[x]]);

(* Formula 397: branch-correct.  The classical form 2 (Sin[x/2] + Cos[x/2])
   differentiates to Cos[x/2] - Sin[x/2], matching Sqrt[1 - Sin[x]] only
   when Cos[x/2] - Sin[x/2] >= 0.  Keeping the radical literally,
   2 Cos[x]/Sqrt[1-Sin[x]] differentiates to (Sin^2 - 2 Sin + 1)/(1-Sin)^(3/2)
   = Sqrt[1-Sin[x]] on every branch. *)
IntegrateTable[Sqrt[1 - Sin[x_]], x_] := (2 Cos[x])/Sqrt[1 - Sin[x]];

(* Formulas 398-401: branch-correct.  The classical primary-branch forms
   Sqrt[2] Log[Tan[x/4 (+ shift)]] are real on only every other inter-pole
   interval of the integrand -- on the rest Tan[.] turns negative and the Log
   goes complex (and even where real its sign is wrong, because the integrand
   carries an absolute value, e.g. 1/Sqrt[1-Cos[x]] = 1/(Sqrt[2] |Sin[x/2]|)).
   These integrands have poles (1/Sqrt[1-Cos] at x=2k Pi, etc.), so no
   continuous antiderivative spans a pole; branch-correctness means the
   derivative equals the integrand on every pole-bounded interval.  The forms
   below keep the integrand's radical literally: the factor
   Sqrt[1-Cos[x]]/(Sqrt[2] Sin[x/2]) = sign(Sin[x/2]) supplies the sign the
   ArcTanh primitive of the signed reciprocal otherwise drops, so the
   derivative matches on every branch (incl. negative x). *)

(* Formula 398 *)
IntegrateTable[1/Sqrt[1 - Cos[x_]], x_] := -ArcTanh[Cos[x/2]] Csc[x/2] Sqrt[1 - Cos[x]];

(* Formula 399 *)
IntegrateTable[1/Sqrt[1 + Cos[x_]], x_] := ArcTanh[Sin[x/2]] Sec[x/2] Sqrt[1 + Cos[x]];

(* Formula 400 *)
IntegrateTable[1/Sqrt[1 - Sin[x_]], x_] := -ArcTanh[Cos[x/2 - Pi/4]] Csc[x/2 - Pi/4] Sqrt[1 - Sin[x]];

(* Formula 401 *)
IntegrateTable[1/Sqrt[1 + Sin[x_]], x_] := ArcTanh[Sin[x/2 - Pi/4]] Sec[x/2 - Pi/4] Sqrt[1 + Sin[x]];

(* Formula 402 *)
IntegrateTable[Tan[a_. x_]^2, x_] /; FreeQ[a, x] := 
  Tan[a x]/a - x;

(* Formula 403 *)
IntegrateTable[Tan[a_. x_]^3, x_] /; FreeQ[a, x] := 
  Tan[a x]^2/(2 a) + 1/a Log[Cos[a x]];

(* Formula 404 *)
IntegrateTable[Tan[a_. x_]^4, x_] /; FreeQ[a, x] := 
  Tan[a x]^3/(3 a) - Tan[a x]/a + x;

(* Formula 405 *)
IntegrateTable[Tan[a_. x_]^n_, x_] /; FreeQ[{a, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 :=
  Tan[a x]^(n - 1)/(a (n - 1)) - IntegrateTable[Tan[a x]^(n - 2), x];

(* Formula 406 *)
IntegrateTable[Cot[a_. x_]^2, x_] /; FreeQ[a, x] := 
  -Cot[a x]/a - x;

(* Formula 407 *)
IntegrateTable[Cot[a_. x_]^3, x_] /; FreeQ[a, x] := 
  -Cot[a x]^2/(2 a) - 1/a Log[Sin[a x]];

(* Formula 408 *)
IntegrateTable[Cot[a_. x_]^4, x_] /; FreeQ[a, x] := 
  -Cot[a x]^3/(3 a) + Cot[a x]/a + x;

(* Formula 409 *)
IntegrateTable[Cot[a_. x_]^n_, x_] /; FreeQ[{a, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 :=
  -Cot[a x]^(n - 1)/(a (n - 1)) - IntegrateTable[Cot[a x]^(n - 2), x];

(* Formula 410 *)
IntegrateTable[x_/Sin[a_. x_]^2, x_] /; FreeQ[a, x] := 
  -(x Cot[a x])/a + 1/a^2 Log[Sin[a x]];

(* Formula 411 *)
IntegrateTable[x_/Sin[a_. x_]^n_, x_] /; FreeQ[{a, n}, x] && n =!= 1 && n =!= 2 && IntegerQ[n] && n > 2 :=
  -(x Cos[a x])/(a (n - 1) Sin[a x]^(n - 1)) - 1/(a^2 (n - 1) (n - 2) Sin[a x]^(n - 2)) + (n - 2)/(n - 1) IntegrateTable[x/Sin[a x]^(n - 2), x];

(* Formula 412 *)
IntegrateTable[x_/Cos[a_. x_]^2, x_] /; FreeQ[a, x] := 
  (x Tan[a x])/a + 1/a^2 Log[Cos[a x]];

(* Formula 413 *)
IntegrateTable[x_/Cos[a_. x_]^n_, x_] /; FreeQ[{a, n}, x] && n =!= 1 && n =!= 2 && IntegerQ[n] && n > 2 :=
  (x Sin[a x])/(a (n - 1) Cos[a x]^(n - 1)) - 1/(a^2 (n - 1) (n - 2) Cos[a x]^(n - 2)) + (n - 2)/(n - 1) IntegrateTable[x/Cos[a x]^(n - 2), x];

(* Formulas 410-413 keyed on the Csc / Sec HEADS: Mathilda folds x/Sin^n into
   x Csc^n, so the reciprocal forms above never match.  Even n closes on the
   x Csc^2 / x Sec^2 bases; odd n bottoms out on the non-elementary x Csc / x Sec
   (a dilogarithm) and is left unresolved, as with x/Sin, x/Cos. *)
IntegrateTable[x_ Csc[a_. x_]^2, x_] /; FreeQ[a, x] :=
  -x Cot[a x]/a + Log[Sin[a x]]/a^2;
IntegrateTable[x_ Csc[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 2 :=
  -x Cot[a x] Csc[a x]^(n - 2)/(a (n - 1)) - Csc[a x]^(n - 2)/(a^2 (n - 1) (n - 2)) + (n - 2)/(n - 1) IntegrateTable[x Csc[a x]^(n - 2), x];
IntegrateTable[x_ Sec[a_. x_]^2, x_] /; FreeQ[a, x] :=
  x Tan[a x]/a + Log[Cos[a x]]/a^2;
IntegrateTable[x_ Sec[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 2 :=
  x Tan[a x] Sec[a x]^(n - 2)/(a (n - 1)) - Sec[a x]^(n - 2)/(a^2 (n - 1) (n - 2)) + (n - 2)/(n - 1) IntegrateTable[x Sec[a x]^(n - 2), x];

(* Formulas 414-421.  The squared coefficient b_^2 (which never matched a numeric
   argument) is bound linearly as b_ and recovered via Sqrt[b]; a_ is the
   frequency (already linear).  The 1 - b^2 Sin^2 members are matched by the same
   1 + b_ Sin^2 form under TrueQ[b < 0] (Sqrt[-b]), since a `1 - b_ Sin^2` pattern
   will not bind b_ to a negative literal. *)

(* Formula 414 *)
IntegrateTable[Sin[a_. x_]/Sqrt[1 + b_ Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] && Not[TrueQ[b < 0]] :=
  -1/(a Sqrt[b]) ArcSin[(Sqrt[b] Cos[a x])/Sqrt[1 + b]];

(* Formula 415 *)
IntegrateTable[Sin[a_. x_]/Sqrt[1 + b_ Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] && TrueQ[b < 0] :=
  -1/(a Sqrt[-b]) Log[Sqrt[-b] Cos[a x] + Sqrt[1 + b Sin[a x]^2]];

(* Formula 416 *)
IntegrateTable[Sin[a_. x_] Sqrt[1 + b_ Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] && Not[TrueQ[b < 0]] :=
  -(Cos[a x]/(2 a)) Sqrt[1 + b Sin[a x]^2] - (1 + b)/(2 a Sqrt[b]) ArcSin[(Sqrt[b] Cos[a x])/Sqrt[1 + b]];

(* Formula 417 *)
IntegrateTable[Sin[a_. x_] Sqrt[1 + b_ Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] && TrueQ[b < 0] :=
  -(Cos[a x]/(2 a)) Sqrt[1 + b Sin[a x]^2] - (1 + b)/(2 a Sqrt[-b]) Log[Sqrt[-b] Cos[a x] + Sqrt[1 + b Sin[a x]^2]];

(* Formula 418 *)
IntegrateTable[Cos[a_. x_]/Sqrt[1 + b_ Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] && Not[TrueQ[b < 0]] :=
  1/(a Sqrt[b]) Log[Sqrt[b] Sin[a x] + Sqrt[1 + b Sin[a x]^2]];

(* Formula 419 *)
IntegrateTable[Cos[a_. x_]/Sqrt[1 + b_ Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] && TrueQ[b < 0] :=
  1/(a Sqrt[-b]) ArcSin[Sqrt[-b] Sin[a x]];

(* Formula 420 *)
IntegrateTable[Cos[a_. x_] Sqrt[1 + b_ Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] && Not[TrueQ[b < 0]] :=
  (Sin[a x]/(2 a)) Sqrt[1 + b Sin[a x]^2] + 1/(2 a Sqrt[b]) Log[Sqrt[b] Sin[a x] + Sqrt[1 + b Sin[a x]^2]];

(* Formula 421 *)
IntegrateTable[Cos[a_. x_] Sqrt[1 + b_ Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] && TrueQ[b < 0] :=
  (Sin[a x]/(2 a)) Sqrt[1 + b Sin[a x]^2] + 1/(2 a Sqrt[-b]) ArcSin[Sqrt[-b] Sin[a x]];

(* Formula 422: Primary Branch *)
IntegrateTable[1/Sqrt[a_ + b_. Tan[c_. x_]^2], x_] /; FreeQ[{a, b, c}, x] && a > Abs[b] := 
  1/(c Sqrt[a - b]) ArcSin[Sqrt[(a - b)/a] Sin[c x]];

(* Formula 427 *)
IntegrateTable[ArcSin[a_. x_], x_] /; FreeQ[a, x] := 
  x ArcSin[a x] + Sqrt[1 - a^2 x^2]/a;

(* Formula 428 *)
IntegrateTable[ArcCos[a_. x_], x_] /; FreeQ[a, x] := 
  x ArcCos[a x] - Sqrt[1 - a^2 x^2]/a;

(* Formula 429 *)
IntegrateTable[ArcTan[a_. x_], x_] /; FreeQ[a, x] := 
  x ArcTan[a x] - 1/(2 a) Log[1 + a^2 x^2];

(* Formula 430 *)
IntegrateTable[ArcCot[a_. x_], x_] /; FreeQ[a, x] := 
  x ArcCot[a x] + 1/(2 a) Log[1 + a^2 x^2];

(* Formula 431 *)
IntegrateTable[ArcSec[a_. x_], x_] /; FreeQ[a, x] := 
  x ArcSec[a x] - 1/a Log[a x + Sqrt[a^2 x^2 - 1]];

(* Formula 432 *)
IntegrateTable[ArcCsc[a_. x_], x_] /; FreeQ[a, x] := 
  x ArcCsc[a x] + 1/a Log[a x + Sqrt[a^2 x^2 - 1]];

(* Formula 437 *)
IntegrateTable[x_ ArcSin[a_. x_], x_] /; FreeQ[a, x] := 
  1/(4 a^2) ((2 a^2 x^2 - 1) ArcSin[a x] + a x Sqrt[1 - a^2 x^2]);

(* Formula 438 *)
IntegrateTable[x_ ArcCos[a_. x_], x_] /; FreeQ[a, x] := 
  1/(4 a^2) ((2 a^2 x^2 - 1) ArcCos[a x] - a x Sqrt[1 - a^2 x^2]);

(* Formula 439 *)
IntegrateTable[x_^n_ ArcSin[a_. x_], x_] /; FreeQ[{a, n}, x] && n =!= -1 := 
  x^(n + 1)/(n + 1) ArcSin[a x] - a/(n + 1) IntegrateTable[x^(n + 1)/Sqrt[1 - a^2 x^2], x];

(* Formula 440 *)
IntegrateTable[x_^n_ ArcCos[a_. x_], x_] /; FreeQ[{a, n}, x] && n =!= -1 := 
  x^(n + 1)/(n + 1) ArcCos[a x] + a/(n + 1) IntegrateTable[x^(n + 1)/Sqrt[1 - a^2 x^2], x];

(* Formula 441 *)
IntegrateTable[x_ ArcTan[a_. x_], x_] /; FreeQ[a, x] := 
  (1 + a^2 x^2)/(2 a^2) ArcTan[a x] - x/(2 a);

(* Formula 442 *)
IntegrateTable[x_^n_ ArcTan[a_. x_], x_] /; FreeQ[{a, n}, x] && n =!= -1 := 
  x^(n + 1)/(n + 1) ArcTan[a x] - a/(n + 1) IntegrateTable[x^(n + 1)/(1 + a^2 x^2), x];

(* Formula 443 *)
IntegrateTable[x_ ArcCot[a_. x_], x_] /; FreeQ[a, x] := 
  (1 + a^2 x^2)/(2 a^2) ArcCot[a x] + x/(2 a);

(* Formula 444 *)
IntegrateTable[x_^n_ ArcCot[a_. x_], x_] /; FreeQ[{a, n}, x] && n =!= -1 := 
  x^(n + 1)/(n + 1) ArcCot[a x] + a/(n + 1) IntegrateTable[x^(n + 1)/(1 + a^2 x^2), x];

(* Formula 445 *)
IntegrateTable[ArcSin[a_. x_]/x_^2, x_] /; FreeQ[a, x] := 
  a Log[(1 - Sqrt[1 - a^2 x^2])/x] - ArcSin[a x]/x;

(* Formula 446 *)
IntegrateTable[ArcCos[a_. x_]/x_^2, x_] /; FreeQ[a, x] := 
  -ArcCos[a x]/x + a Log[(1 + Sqrt[1 - a^2 x^2])/x];

(* Formula 447 *)
IntegrateTable[ArcTan[a_. x_]/x_^2, x_] /; FreeQ[a, x] := 
  -ArcTan[a x]/x - a/2 Log[(1 + a^2 x^2)/x^2];

(* Formula 448 *)
IntegrateTable[ArcCot[a_. x_]/x_^2, x_] /; FreeQ[a, x] := 
  -ArcCot[a x]/x - a/2 Log[x^2/(1 + a^2 x^2)];

(* Formula 449 *)
IntegrateTable[ArcSin[a_. x_]^2, x_] /; FreeQ[a, x] := 
  x ArcSin[a x]^2 - 2 x + (2 Sqrt[1 - a^2 x^2])/a ArcSin[a x];

(* Formula 450 *)
IntegrateTable[ArcCos[a_. x_]^2, x_] /; FreeQ[a, x] := 
  x ArcCos[a x]^2 - 2 x - (2 Sqrt[1 - a^2 x^2])/a ArcCos[a x];

(* Formula 453.  The quadratic coefficient is matched as a_ (= -b^2) and the *)
(* linear coefficient is recovered as Sqrt[-a] on the RHS: matching a_^2 in    *)
(* the pattern would require the matcher to invert a square, which it does not *)
(* do, so the coefficient is bound linearly and re-squared via the Condition.  *)
IntegrateTable[ArcSin[b_. x_]/Sqrt[c_ + a_. x_^2], x_] /; FreeQ[{a, b, c}, x] && c === 1 && a === -b^2 :=
  ArcSin[Sqrt[-a] x]^2/(2 Sqrt[-a]);

(* Formula 454: Mapped explicitly from the integral recurrence relation.  The  *)
(* optional exponent x_^n_. also matches the bare x (n = 1), where the trailing *)
(* recurrence term vanishes (coefficient n - 1 = 0), so odd powers bottom out   *)
(* cleanly without a separate n = 1 rule. *)
IntegrateTable[x_^n_. ArcSin[b_. x_]/Sqrt[c_ + a_. x_^2], x_] /; FreeQ[{a, b, c, n}, x] && c === 1 && a === -b^2 && IntegerQ[n] && n > 0 :=
  -(x^(n - 1)/(n (-a))) Sqrt[1 + a x^2] ArcSin[Sqrt[-a] x] + x^n/(n^2 Sqrt[-a]) + (n - 1)/(n (-a)) IntegrateTable[(x^(n - 2) ArcSin[Sqrt[-a] x])/Sqrt[1 + a x^2], x];

(* Formula 455 *)
IntegrateTable[ArcCos[b_. x_]/Sqrt[c_ + a_. x_^2], x_] /; FreeQ[{a, b, c}, x] && c === 1 && a === -b^2 :=
  -ArcCos[Sqrt[-a] x]^2/(2 Sqrt[-a]);

(* Formula 456 *)
IntegrateTable[x_^n_. ArcCos[b_. x_]/Sqrt[c_ + a_. x_^2], x_] /; FreeQ[{a, b, c, n}, x] && c === 1 && a === -b^2 && IntegerQ[n] && n > 0 :=
  -(x^(n - 1)/(n (-a))) Sqrt[1 + a x^2] ArcCos[Sqrt[-a] x] - x^n/(n^2 Sqrt[-a]) + (n - 1)/(n (-a)) IntegrateTable[(x^(n - 2) ArcCos[Sqrt[-a] x])/Sqrt[1 + a x^2], x];

(* Formula 457 *)
IntegrateTable[ArcTan[b_. x_]/(c_ + a_. x_^2), x_] /; FreeQ[{a, b, c}, x] && c === 1 && a === b^2 :=
  ArcTan[Sqrt[a] x]^2/(2 Sqrt[a]);

(* Formula 458 *)
IntegrateTable[ArcCot[b_. x_]/(c_ + a_. x_^2), x_] /; FreeQ[{a, b, c}, x] && c === 1 && a === b^2 :=
  -ArcCot[Sqrt[a] x]^2/(2 Sqrt[a]);

(* Formula 459 *)
IntegrateTable[x_ ArcSec[a_. x_], x_] /; FreeQ[a, x] := 
  x^2/2 ArcSec[a x] - 1/(2 a^2) Sqrt[a^2 x^2 - 1];

(* Formula 460 *)
IntegrateTable[x_^n_ ArcSec[a_. x_], x_] /; FreeQ[{a, n}, x] && n =!= -1 := 
  x^(n + 1)/(n + 1) ArcSec[a x] - 1/(n + 1) IntegrateTable[x^n/Sqrt[a^2 x^2 - 1], x];

(* Formula 461 *)
IntegrateTable[ArcSec[a_. x_]/x_^2, x_] /; FreeQ[a, x] := 
  -ArcSec[a x]/x + Sqrt[a^2 x^2 - 1]/x;

(* Formula 462 *)
IntegrateTable[x_ ArcCsc[a_. x_], x_] /; FreeQ[a, x] := 
  x^2/2 ArcCsc[a x] + 1/(2 a^2) Sqrt[a^2 x^2 - 1];

(* Formula 463 *)
IntegrateTable[x_^n_ ArcCsc[a_. x_], x_] /; FreeQ[{a, n}, x] && n =!= -1 := 
  x^(n + 1)/(n + 1) ArcCsc[a x] + 1/(n + 1) IntegrateTable[x^n/Sqrt[a^2 x^2 - 1], x];

(* Formula 464 *)
IntegrateTable[ArcCsc[a_. x_]/x_^2, x_] /; FreeQ[a, x] :=
  -ArcCsc[a x]/x - Sqrt[a^2 x^2 - 1]/x;

(* ============================================================ *)
(* Inverse-hyperbolic analogs of Formulas 427-464.              *)
(* Each rule below mirrors the inverse-trig entry immediately   *)
(* above, obtained by the substitutions                         *)
(*   ArcSin -> ArcSinh, ArcCos -> ArcCosh, ArcTan -> ArcTanh,   *)
(*   ArcCot -> ArcCoth, ArcSec -> ArcSech, ArcCsc -> ArcCsch,   *)
(* with the radicand / denominator sign flips forced by the     *)
(* inverse-hyperbolic derivatives:                              *)
(*   D[ArcSinh[u]] = 1/Sqrt[1+u^2],  D[ArcCosh[u]] = 1/Sqrt[u^2-1], *)
(*   D[ArcTanh[u]] = D[ArcCoth[u]] = 1/(1-u^2).                  *)
(* Every rule is diff-back verified (D[F,x]-integrand -> 0) by  *)
(* the CRC corpus runner; the ArcSech/ArcCsch entries carry a   *)
(* |x| branch factor that closes numerically on the positive    *)
(* real axis.                                                   *)
(* ============================================================ *)

(* Formula 427h *)
IntegrateTable[ArcSinh[a_. x_], x_] /; FreeQ[a, x] :=
  x ArcSinh[a x] - Sqrt[1 + a^2 x^2]/a;

(* Formula 428h *)
IntegrateTable[ArcCosh[a_. x_], x_] /; FreeQ[a, x] :=
  x ArcCosh[a x] - Sqrt[a^2 x^2 - 1]/a;

(* Formula 429h *)
IntegrateTable[ArcTanh[a_. x_], x_] /; FreeQ[a, x] :=
  x ArcTanh[a x] + 1/(2 a) Log[1 - a^2 x^2];

(* Formula 430h *)
IntegrateTable[ArcCoth[a_. x_], x_] /; FreeQ[a, x] :=
  x ArcCoth[a x] + 1/(2 a) Log[a^2 x^2 - 1];

(* Formula 431h *)
IntegrateTable[ArcSech[a_. x_], x_] /; FreeQ[a, x] :=
  x ArcSech[a x] + 1/a ArcSin[a x];

(* Formula 432h *)
IntegrateTable[ArcCsch[a_. x_], x_] /; FreeQ[a, x] :=
  x ArcCsch[a x] + 1/a ArcSinh[a x];

(* Formula 437h *)
IntegrateTable[x_ ArcSinh[a_. x_], x_] /; FreeQ[a, x] :=
  1/(4 a^2) ((2 a^2 x^2 + 1) ArcSinh[a x] - a x Sqrt[1 + a^2 x^2]);

(* Formula 438h *)
IntegrateTable[x_ ArcCosh[a_. x_], x_] /; FreeQ[a, x] :=
  1/(4 a^2) ((2 a^2 x^2 - 1) ArcCosh[a x] - a x Sqrt[a^2 x^2 - 1]);

(* Formula 439h *)
IntegrateTable[x_^n_ ArcSinh[a_. x_], x_] /; FreeQ[{a, n}, x] && n =!= -1 :=
  x^(n + 1)/(n + 1) ArcSinh[a x] - a/(n + 1) IntegrateTable[x^(n + 1)/Sqrt[1 + a^2 x^2], x];

(* Formula 440h *)
IntegrateTable[x_^n_ ArcCosh[a_. x_], x_] /; FreeQ[{a, n}, x] && n =!= -1 :=
  x^(n + 1)/(n + 1) ArcCosh[a x] - a/(n + 1) IntegrateTable[x^(n + 1)/Sqrt[a^2 x^2 - 1], x];

(* Formula 441h *)
IntegrateTable[x_ ArcTanh[a_. x_], x_] /; FreeQ[a, x] :=
  (a^2 x^2 - 1)/(2 a^2) ArcTanh[a x] + x/(2 a);

(* Formula 442h *)
IntegrateTable[x_^n_ ArcTanh[a_. x_], x_] /; FreeQ[{a, n}, x] && n =!= -1 :=
  x^(n + 1)/(n + 1) ArcTanh[a x] - a/(n + 1) IntegrateTable[x^(n + 1)/(1 - a^2 x^2), x];

(* Formula 443h *)
IntegrateTable[x_ ArcCoth[a_. x_], x_] /; FreeQ[a, x] :=
  (a^2 x^2 - 1)/(2 a^2) ArcCoth[a x] + x/(2 a);

(* Formula 444h *)
IntegrateTable[x_^n_ ArcCoth[a_. x_], x_] /; FreeQ[{a, n}, x] && n =!= -1 :=
  x^(n + 1)/(n + 1) ArcCoth[a x] - a/(n + 1) IntegrateTable[x^(n + 1)/(1 - a^2 x^2), x];

(* Formula 445h *)
IntegrateTable[ArcSinh[a_. x_]/x_^2, x_] /; FreeQ[a, x] :=
  -ArcSinh[a x]/x - a Log[(1 + Sqrt[1 + a^2 x^2])/x];

(* Formula 446h *)
IntegrateTable[ArcCosh[a_. x_]/x_^2, x_] /; FreeQ[a, x] :=
  -ArcCosh[a x]/x + a ArcTan[Sqrt[a^2 x^2 - 1]];

(* Formula 447h *)
IntegrateTable[ArcTanh[a_. x_]/x_^2, x_] /; FreeQ[a, x] :=
  -ArcTanh[a x]/x - a/2 Log[(1 - a^2 x^2)/x^2];

(* Formula 448h *)
IntegrateTable[ArcCoth[a_. x_]/x_^2, x_] /; FreeQ[a, x] :=
  -ArcCoth[a x]/x + a/2 Log[x^2/(1 - a^2 x^2)];

(* Formula 449h *)
IntegrateTable[ArcSinh[a_. x_]^2, x_] /; FreeQ[a, x] :=
  x ArcSinh[a x]^2 + 2 x - (2 Sqrt[1 + a^2 x^2])/a ArcSinh[a x];

(* Formula 450h *)
IntegrateTable[ArcCosh[a_. x_]^2, x_] /; FreeQ[a, x] :=
  x ArcCosh[a x]^2 + 2 x - (2 Sqrt[a^2 x^2 - 1])/a ArcCosh[a x];

(* Formula 453h *)
IntegrateTable[ArcSinh[b_. x_]/Sqrt[c_ + a_. x_^2], x_] /; FreeQ[{a, b, c}, x] && c === 1 && a === b^2 :=
  ArcSinh[Sqrt[a] x]^2/(2 Sqrt[a]);

(* Formula 454h *)
IntegrateTable[x_^n_. ArcSinh[b_. x_]/Sqrt[c_ + a_. x_^2], x_] /; FreeQ[{a, b, c, n}, x] && c === 1 && a === b^2 && IntegerQ[n] && n > 0 :=
  (x^(n - 1)/(n a)) Sqrt[1 + a x^2] ArcSinh[Sqrt[a] x] - x^n/(n^2 Sqrt[a]) - (n - 1)/(n a) IntegrateTable[(x^(n - 2) ArcSinh[Sqrt[a] x])/Sqrt[1 + a x^2], x];

(* Formula 455h *)
IntegrateTable[ArcCosh[b_. x_]/Sqrt[c_ + a_. x_^2], x_] /; FreeQ[{a, b, c}, x] && c === -1 && a === b^2 :=
  ArcCosh[Sqrt[a] x]^2/(2 Sqrt[a]);

(* Formula 456h *)
IntegrateTable[x_^n_. ArcCosh[b_. x_]/Sqrt[c_ + a_. x_^2], x_] /; FreeQ[{a, b, c, n}, x] && c === -1 && a === b^2 && IntegerQ[n] && n > 0 :=
  (x^(n - 1)/(n a)) Sqrt[a x^2 - 1] ArcCosh[Sqrt[a] x] - x^n/(n^2 Sqrt[a]) + (n - 1)/(n a) IntegrateTable[(x^(n - 2) ArcCosh[Sqrt[a] x])/Sqrt[a x^2 - 1], x];

(* Formula 457h *)
IntegrateTable[ArcTanh[b_. x_]/(c_ + a_. x_^2), x_] /; FreeQ[{a, b, c}, x] && c === 1 && a === -b^2 :=
  ArcTanh[Sqrt[-a] x]^2/(2 Sqrt[-a]);

(* Formula 458h *)
IntegrateTable[ArcCoth[b_. x_]/(c_ + a_. x_^2), x_] /; FreeQ[{a, b, c}, x] && c === 1 && a === -b^2 :=
  ArcCoth[Sqrt[-a] x]^2/(2 Sqrt[-a]);

(* Formula 459h *)
IntegrateTable[x_ ArcSech[a_. x_], x_] /; FreeQ[a, x] :=
  x^2/2 ArcSech[a x] - 1/(2 a^2) Sqrt[1 - a^2 x^2];

(* Formula 460h *)
IntegrateTable[x_^n_ ArcSech[a_. x_], x_] /; FreeQ[{a, n}, x] && n =!= -1 :=
  x^(n + 1)/(n + 1) ArcSech[a x] + 1/(n + 1) IntegrateTable[x^n/Sqrt[1 - a^2 x^2], x];

(* Formula 461h *)
IntegrateTable[ArcSech[a_. x_]/x_^2, x_] /; FreeQ[a, x] :=
  -ArcSech[a x]/x + Sqrt[1 - a^2 x^2]/x;

(* Formula 462h *)
IntegrateTable[x_ ArcCsch[a_. x_], x_] /; FreeQ[a, x] :=
  x^2/2 ArcCsch[a x] + 1/(2 a^2) Sqrt[1 + a^2 x^2];

(* Formula 463h *)
IntegrateTable[x_^n_ ArcCsch[a_. x_], x_] /; FreeQ[{a, n}, x] && n =!= -1 :=
  x^(n + 1)/(n + 1) ArcCsch[a x] + 1/(n + 1) IntegrateTable[x^n/Sqrt[1 + a^2 x^2], x];

(* Formula 464h *)
IntegrateTable[ArcCsch[a_. x_]/x_^2, x_] /; FreeQ[a, x] :=
  -ArcCsch[a x]/x + Sqrt[1 + a^2 x^2]/x;

(* Formula 465 *)
IntegrateTable[Log[x_], x_] := 
  x Log[x] - x;

(* Formula 466 *)
IntegrateTable[x_ Log[x_], x_] := 
  x^2/2 Log[x] - x^2/4;

(* Formula 467 *)
IntegrateTable[x_^2 Log[x_], x_] := 
  x^3/3 Log[x] - x^3/9;

(* Formula 468 *)
IntegrateTable[x_^n_ Log[x_], x_] /; FreeQ[n, x] && n =!= -1 := 
  x^(n + 1)/(n + 1) Log[x] - x^(n + 1)/(n + 1)^2;

(* Formula 469 *)
IntegrateTable[Log[x_]^2, x_] := 
  x Log[x]^2 - 2 x Log[x] + 2 x;

(* Formula 471 *)
IntegrateTable[Log[x_]^n_/x_, x_] /; FreeQ[n, x] && n =!= -1 := 
  Log[x]^(n + 1)/(n + 1);

(* Formula 473 *)
IntegrateTable[1/(x_ Log[x_]), x_] := 
  Log[Log[x]];

(* Formula 474 *)
IntegrateTable[1/(x_ Log[x_]^n_), x_] /; FreeQ[n, x] && n =!= 1 := 
  1/((1 - n) Log[x]^(n - 1));

(* Formula 475 *)
IntegrateTable[x_^m_/Log[x_]^n_, x_] /; FreeQ[{m, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 :=
  x^(m + 1)/((1 - n) Log[x]^(n - 1)) + (m + 1)/(n - 1) IntegrateTable[x^m/Log[x]^(n - 1), x];

(* Formula 476 *)
IntegrateTable[x_^m_ Log[x_]^n_, x_] /; FreeQ[{m, n}, x] && IntegerQ[n] && n > 0 :=
  (x^(m + 1) Log[x]^n)/(m + 1) - n/(m + 1) IntegrateTable[x^m Log[x]^(n - 1), x];

(* Formula 477 *)
IntegrateTable[x_^p_ Cos[b_. Log[x_]], x_] /; FreeQ[{p, b}, x] := 
  x^(p + 1)/((p + 1)^2 + b^2) (b Sin[b Log[x]] + (p + 1) Cos[b Log[x]]);

(* Formula 478 *)
IntegrateTable[x_^p_ Sin[b_. Log[x_]], x_] /; FreeQ[{p, b}, x] := 
  x^(p + 1)/((p + 1)^2 + b^2) ((p + 1) Sin[b Log[x]] - b Cos[b Log[x]]);

(* Formula 479 *)
IntegrateTable[Log[a_. x_ + b_], x_] /; FreeQ[{a, b}, x] := 
  (a x + b)/a Log[a x + b] - x;

(* Formula 480 *)
IntegrateTable[Log[a_. x_ + b_]/x_^2, x_] /; FreeQ[{a, b}, x] := 
  a/b Log[x] - (a x + b)/(b x) Log[a x + b];

(* Formula 483 *)
IntegrateTable[Log[(x_ + a_)/(x_ - a_)], x_] /; FreeQ[a, x] := 
  (x + a) Log[x + a] - (x - a) Log[x - a];

(* Formula 485 *)
IntegrateTable[1/x_^2 Log[(x_ + a_)/(x_ - a_)], x_] /; FreeQ[a, x] := 
  1/x Log[(x - a)/(x + a)] - 1/a Log[(x^2 - a^2)/x^2];

(* Formula 486 (Mapped using With for X variables) *)
IntegrateTable[Log[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c}, x] && 4 a c - b^2 > 0 := 
  With[{q = 4 a c - b^2, X = a + b x + c x^2}, 
    (x + b/(2 c)) Log[X] - 2 x + Sqrt[q]/c ArcTan[(2 c x + b)/Sqrt[q]]
  ];
IntegrateTable[Log[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c}, x] && b^2 - 4 a c > 0 := 
  With[{q = b^2 - 4 a c, X = a + b x + c x^2}, 
    (x + b/(2 c)) Log[X] - 2 x + Sqrt[q]/c ArcTanh[(2 c x + b)/Sqrt[q]]
  ];

(* Formula 487 *)
IntegrateTable[x_^n_ Log[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c, n}, x] && n =!= -1 := 
  With[{X = a + b x + c x^2}, 
    x^(n + 1)/(n + 1) Log[X] - (2 c)/(n + 1) IntegrateTable[x^(n + 2)/X, x] - b/(n + 1) IntegrateTable[x^(n + 1)/X, x]
  ];

(* Formulas 488-495.  These rules matched a squared constant as `a_^2`,
   which the pattern matcher cannot bind against a numeric argument (it
   does not invert the square), so they never fired for concrete inputs.
   Rewritten to bind the constant linearly as `a_` and re-express the
   linear part via Sqrt[a].  Because `x_^2 - a_` will not bind `a_` to a
   negative literal, BOTH signs are matched with the single `x_^2 + a_`
   form, disambiguated by the guard: `Not[TrueQ[a < 0]]` for the plus rule
   (also admits a symbolic `a`), `TrueQ[a < 0]` for the minus rule (using
   Sqrt[-a]).  When the RHS is even in the original constant (only a^2
   appears, e.g. 490-493) the negative-a plus form already reproduces the
   minus case, so a single unguarded rule covers both. *)

(* Formula 488 *)
IntegrateTable[Log[x_^2 + a_], x_] /; FreeQ[a, x] && Not[TrueQ[a < 0]] :=
  x Log[x^2 + a] - 2 x + 2 Sqrt[a] ArcTan[x/Sqrt[a]];

(* Formula 489 *)
IntegrateTable[Log[x_^2 + a_], x_] /; FreeQ[a, x] && TrueQ[a < 0] :=
  x Log[x^2 + a] - 2 x + Sqrt[-a] Log[(x + Sqrt[-a])/(x - Sqrt[-a])];

(* Formula 490 *)
IntegrateTable[x_ Log[x_^2 + a_], x_] /; FreeQ[a, x] :=
  1/2 (x^2 + a) Log[x^2 + a] - x^2/2;

(* Formula 491 *)
IntegrateTable[Log[x_ + Sqrt[x_^2 + a_]], x_] /; FreeQ[a, x] :=
  x Log[x + Sqrt[x^2 + a]] - Sqrt[x^2 + a];

(* Formula 492 *)
IntegrateTable[x_ Log[x_ + Sqrt[x_^2 + a_]], x_] /; FreeQ[a, x] :=
  (x^2/2 + a/4) Log[x + Sqrt[x^2 + a]] - (x Sqrt[x^2 + a])/4;

(* Formula 493 *)
IntegrateTable[x_^m_ Log[x_ + Sqrt[x_^2 + a_]], x_] /; FreeQ[{a, m}, x] :=
  x^(m + 1)/(m + 1) Log[x + Sqrt[x^2 + a]] - 1/(m + 1) IntegrateTable[x^(m + 1)/Sqrt[x^2 + a], x];

(* Formula 494 *)
IntegrateTable[Log[x_ + Sqrt[x_^2 + a_]]/x_^2, x_] /; FreeQ[a, x] && Not[TrueQ[a < 0]] :=
  -Log[x + Sqrt[x^2 + a]]/x - 1/Sqrt[a] Log[(Sqrt[a] + Sqrt[x^2 + a])/x];

(* Formula 495 *)
IntegrateTable[Log[x_ + Sqrt[x_^2 + a_]]/x_^2, x_] /; FreeQ[a, x] && TrueQ[a < 0] :=
  -Log[x + Sqrt[x^2 + a]]/x + 1/Sqrt[-a] ArcSec[x/Sqrt[-a]];

(* Formula 497 *)
IntegrateTable[E^x_, x_] := E^x;

(* Formula 498 *)
IntegrateTable[E^(-x_), x_] := -E^(-x);

(* Formula 499 *)
IntegrateTable[E^(a_. x_), x_] /; FreeQ[a, x] := E^(a x)/a;

(* Formula 500 *)
IntegrateTable[x_ E^(a_. x_), x_] /; FreeQ[a, x] := E^(a x)/a^2 (a x - 1);

(* Formula 501: Recursive reduction *)
IntegrateTable[x_^m_ E^(a_. x_), x_] /; FreeQ[{a, m}, x] && m > 0 := 
  (x^m E^(a x))/a - m/a IntegrateTable[x^(m - 1) E^(a x), x];

(* Formula 503 *)
IntegrateTable[E^(a_. x_)/x_^m_, x_] /; FreeQ[{a, m}, x] && m =!= 1 && IntegerQ[m] && m > 1 :=
  1/(1 - m) E^(a x)/x^(m - 1) + a/(m - 1) IntegrateTable[E^(a x)/x^(m - 1), x];

(* Formula 504 *)
IntegrateTable[E^(a_. x_) Log[x_], x_] /; FreeQ[a, x] := 
  (E^(a x) Log[x])/a - 1/a IntegrateTable[E^(a x)/x, x];

(* Formula 505 *)
IntegrateTable[1/(1 + E^x_), x_] := 
  x - Log[1 + E^x];

(* Formula 506 *)
IntegrateTable[1/(a_ + b_. E^(p_. x_)), x_] /; FreeQ[{a, b, p}, x] := 
  x/a - 1/(a p) Log[a + b E^(p x)];

(* Formula 507 *)
IntegrateTable[1/(a_ E^(m_. x_) + b_ E^(-m_. x_)), x_] /; FreeQ[{a, b, m}, x] && a > 0 && b > 0 := 
  1/(m Sqrt[a b]) ArcTan[E^(m x) Sqrt[a/b]];

(* Formula 508 *)
IntegrateTable[1/(a_ E^(m_. x_) - b_ E^(-m_. x_)), x_] /; FreeQ[{a, b, m}, x] && a > 0 && b > 0 := 
  1/(2 m Sqrt[a b]) Log[(Sqrt[a] E^(m x) - Sqrt[b])/(Sqrt[a] E^(m x) + Sqrt[b])];

(* Formula 509 *)
IntegrateTable[a_^x_ - a_^(-x_), x_] /; FreeQ[a, x] && a > 0 := 
  (a^x + a^(-x))/Log[a];

(* Formula 510 *)
IntegrateTable[E^(a_. x_)/(b_ + c_. E^(a_. x_)), x_] /; FreeQ[{a, b, c}, x] := 
  1/(a c) Log[b + c E^(a x)];

(* Formula 511 *)
IntegrateTable[(x_ E^(a_. x_))/(1 + a_. x_)^2, x_] /; FreeQ[a, x] := 
  E^(a x)/(a^2 (1 + a x));

(* Formula 512 *)
IntegrateTable[x_ E^(-x_^2), x_] := 
  -1/2 E^(-x^2);

(* Formula 513 *)
IntegrateTable[E^(a_. x_) Sin[b_. x_], x_] /; FreeQ[{a, b}, x] := 
  E^(a x) (a Sin[b x] - b Cos[b x])/(a^2 + b^2);

(* Formula 514 *)
IntegrateTable[E^(a_. x_) Sin[b_. x_] Sin[c_. x_], x_] /; FreeQ[{a, b, c}, x] := 
  (E^(a x) ((b - c) Sin[(b - c) x] + a Cos[(b - c) x]))/(2 (a^2 + (b - c)^2)) - 
  (E^(a x) ((b + c) Sin[(b + c) x] + a Cos[(b + c) x]))/(2 (a^2 + (b + c)^2));

(* Formula 515 *)
IntegrateTable[E^(a_. x_) Sin[b_. x_] Cos[c_. x_], x_] /; FreeQ[{a, b, c}, x] := 
  (E^(a x) (a Sin[(b - c) x] - (b - c) Cos[(b - c) x]))/(2 (a^2 + (b - c)^2)) + 
  (E^(a x) (a Sin[(b + c) x] - (b + c) Cos[(b + c) x]))/(2 (a^2 + (b + c)^2));

(* Formula 516 *)
IntegrateTable[E^(a_. x_) Sin[b_. x_] Sin[b_. x_ + c_], x_] /; FreeQ[{a, b, c}, x] := 
  (E^(a x) Cos[c])/(2 a) - (E^(a x) (a Cos[2 b x + c] + 2 b Sin[2 b x + c]))/(2 (a^2 + 4 b^2));

(* Formula 517 *)
IntegrateTable[E^(a_. x_) Sin[b_. x_] Cos[b_. x_ + c_], x_] /; FreeQ[{a, b, c}, x] := 
  -(E^(a x) Sin[c])/(2 a) + (E^(a x) (a Sin[2 b x + c] - 2 b Cos[2 b x + c]))/(2 (a^2 + 4 b^2));

(* Formula 518 *)
IntegrateTable[E^(a_. x_) Cos[b_. x_], x_] /; FreeQ[{a, b}, x] := 
  E^(a x)/(a^2 + b^2) (a Cos[b x] + b Sin[b x]);

(* Formula 519 *)
IntegrateTable[E^(a_. x_) Cos[b_. x_] Cos[c_. x_], x_] /; FreeQ[{a, b, c}, x] := 
  (E^(a x) ((b - c) Sin[(b - c) x] + a Cos[(b - c) x]))/(2 (a^2 + (b - c)^2)) + 
  (E^(a x) ((b + c) Sin[(b + c) x] + a Cos[(b + c) x]))/(2 (a^2 + (b + c)^2));

(* Formula 520 *)
IntegrateTable[E^(a_. x_) Cos[b_. x_] Cos[b_. x_ + c_], x_] /; FreeQ[{a, b, c}, x] := 
  (E^(a x) Cos[c])/(2 a) + (E^(a x) (a Cos[2 b x + c] + 2 b Sin[2 b x + c]))/(2 (a^2 + 4 b^2));

(* Formula 521 *)
IntegrateTable[E^(a_. x_) Cos[b_. x_] Sin[b_. x_ + c_], x_] /; FreeQ[{a, b, c}, x] := 
  (E^(a x) Sin[c])/(2 a) + (E^(a x) (a Sin[2 b x + c] - 2 b Cos[2 b x + c]))/(2 (a^2 + 4 b^2));

(* Formula 522 *)
IntegrateTable[E^(a_. x_) Sin[b_. x_]^n_, x_] /; FreeQ[{a, b, n}, x] && n =!= 0 && IntegerQ[n] && n > 1 :=
  1/(a^2 + n^2 b^2) (E^(a x) Sin[b x]^(n - 1) (a Sin[b x] - n b Cos[b x]) + n (n - 1) b^2 IntegrateTable[E^(a x) Sin[b x]^(n - 2), x]);

(* Formula 523 *)
IntegrateTable[E^(a_. x_) Cos[b_. x_]^n_, x_] /; FreeQ[{a, b, n}, x] && n =!= 0 && IntegerQ[n] && n > 1 :=
  1/(a^2 + n^2 b^2) (E^(a x) Cos[b x]^(n - 1) (a Cos[b x] + n b Sin[b x]) + n (n - 1) b^2 IntegrateTable[E^(a x) Cos[b x]^(n - 2), x]);

(* Formula 524 *)
IntegrateTable[x_^m_ E^x_ Sin[x_], x_] /; FreeQ[m, x] && m > 0 := 
  1/2 x^m E^x (Sin[x] - Cos[x]) - m/2 IntegrateTable[x^(m - 1) E^x Sin[x], x] + m/2 IntegrateTable[x^(m - 1) E^x Cos[x], x];

(* Formula 525 *)
IntegrateTable[x_^m_ E^(a_. x_) Sin[b_. x_], x_] /; FreeQ[{a, b, m}, x] && m > 0 := 
  (x^m E^(a x) (a Sin[b x] - b Cos[b x]))/(a^2 + b^2) - m/(a^2 + b^2) IntegrateTable[x^(m - 1) E^(a x) (a Sin[b x] - b Cos[b x]), x];

(* Formula 526 *)
IntegrateTable[x_^m_ E^x_ Cos[x_], x_] /; FreeQ[m, x] && m > 0 := 
  1/2 x^m E^x (Sin[x] + Cos[x]) - m/2 IntegrateTable[x^(m - 1) E^x Sin[x], x] - m/2 IntegrateTable[x^(m - 1) E^x Cos[x], x];

(* Formula 527 *)
IntegrateTable[x_^m_ E^(a_. x_) Cos[b_. x_], x_] /; FreeQ[{a, b, m}, x] && m > 0 := 
  (x^m E^(a x) (a Cos[b x] + b Sin[b x]))/(a^2 + b^2) - m/(a^2 + b^2) IntegrateTable[x^(m - 1) E^(a x) (a Cos[b x] + b Sin[b x]), x];

(* Formula 528 (E^(a x) Cos^m Sin^n) was left as a truncated, RHS-less
   statement in the CRC port and never fired; removed rather than kept as dead
   syntax.  Such products are reached instead through the general cascade. *)

(* E^(a x) times hyperbolic powers.  The circular E^(a x) Sin^n / Cos^n
   reductions above have denominator a^2 + n^2 b^2, which never vanishes for
   real a, b.  The hyperbolic analogue's denominator is a^2 - n^2 b^2, which is
   zero at the resonance a = n b (e.g. E^x Cosh[x], a = b = 1: there e^x =
   cosh x + sinh x forces a secular x/2 term).  The reductions below carry the
   a^2 =!= n^2 b^2 guard and bottom out on the n = 1 bases; the a^2 === b^2
   bases supply the resonant closed form e^(2 a x)/(4 a) +- x/2.  (An exotic
   intermediate resonance a = n b with n > 1, e.g. E^(2x) Cosh[x]^2, is left to
   the general cascade.) *)
IntegrateTable[E^(a_. x_) Sinh[b_. x_], x_] /; FreeQ[{a, b}, x] && a^2 =!= b^2 :=
  E^(a x) (a Sinh[b x] - b Cosh[b x])/(a^2 - b^2);
IntegrateTable[E^(a_. x_) Cosh[b_. x_], x_] /; FreeQ[{a, b}, x] && a^2 =!= b^2 :=
  E^(a x) (a Cosh[b x] - b Sinh[b x])/(a^2 - b^2);
IntegrateTable[E^(a_. x_) Cosh[b_. x_], x_] /; FreeQ[{a, b}, x] && a^2 === b^2 :=
  E^(2 a x)/(4 a) + x/2;
IntegrateTable[E^(a_. x_) Sinh[b_. x_], x_] /; FreeQ[{a, b}, x] && b === a :=
  E^(2 a x)/(4 a) - x/2;
IntegrateTable[E^(a_. x_) Sinh[b_. x_], x_] /; FreeQ[{a, b}, x] && b === -a :=
  -E^(2 a x)/(4 a) + x/2;
IntegrateTable[E^(a_. x_) Sinh[b_. x_]^n_, x_] /; FreeQ[{a, b, n}, x] && IntegerQ[n] && n > 1 && a^2 =!= n^2 b^2 :=
  (E^(a x) Sinh[b x]^(n - 1) (a Sinh[b x] - n b Cosh[b x]))/(a^2 - n^2 b^2) +
  (n (n - 1) b^2)/(a^2 - n^2 b^2) IntegrateTable[E^(a x) Sinh[b x]^(n - 2), x];
IntegrateTable[E^(a_. x_) Cosh[b_. x_]^n_, x_] /; FreeQ[{a, b, n}, x] && IntegerQ[n] && n > 1 && a^2 =!= n^2 b^2 :=
  (E^(a x) Cosh[b x]^(n - 1) (a Cosh[b x] - n b Sinh[b x]))/(a^2 - n^2 b^2) -
  (n (n - 1) b^2)/(a^2 - n^2 b^2) IntegrateTable[E^(a x) Cosh[b x]^(n - 2), x];

(* Formula 534 - 539 *)
(* Single-function bases generalised to argument a x (the CRC originals bound a
   bare x, so every a != 1 hyperbolic reduction that bottomed out here failed to
   close — cf. Sech[2 x]^3). *)
IntegrateTable[Sinh[a_. x_], x_] /; FreeQ[a, x] := Cosh[a x]/a;
IntegrateTable[Cosh[a_. x_], x_] /; FreeQ[a, x] := Sinh[a x]/a;
IntegrateTable[Tanh[a_. x_], x_] /; FreeQ[a, x] := Log[Cosh[a x]]/a;
IntegrateTable[Coth[a_. x_], x_] /; FreeQ[a, x] := Log[Sinh[a x]]/a;
IntegrateTable[Sech[a_. x_], x_] /; FreeQ[a, x] := ArcTan[Sinh[a x]]/a;
IntegrateTable[Csch[a_. x_], x_] /; FreeQ[a, x] := Log[Tanh[(a x)/2]]/a;

(* Hyperbolic sine / cosine integrals (Sinh[a x]/x and Cosh[a x]/x), the
   hyperbolic analogues of the SinIntegral / CosIntegral base cases above.
   D[SinhIntegral[a x], x] = Sinh[a x]/x, D[CoshIntegral[a x], x] = Cosh[a x]/x. *)
IntegrateTable[Sinh[a_. x_]/x_, x_] /; FreeQ[a, x] := SinhIntegral[a x];
IntegrateTable[Cosh[a_. x_]/x_, x_] /; FreeQ[a, x] := CoshIntegral[a x];

(* Formula 540 & 542 (generalised to argument a x) *)
IntegrateTable[x_ Sinh[a_. x_], x_] /; FreeQ[a, x] := x Cosh[a x]/a - Sinh[a x]/a^2;
IntegrateTable[x_ Cosh[a_. x_], x_] /; FreeQ[a, x] := x Sinh[a x]/a - Cosh[a x]/a^2;

(* Formula 541 & 543 (generalised to argument a x) *)
IntegrateTable[x_^n_ Sinh[a_. x_], x_] /; FreeQ[{a, n}, x] && n > 0 :=
  x^n Cosh[a x]/a - n/a IntegrateTable[x^(n - 1) Cosh[a x], x];
IntegrateTable[x_^n_ Cosh[a_. x_], x_] /; FreeQ[{a, n}, x] && n > 0 :=
  x^n Sinh[a x]/a - n/a IntegrateTable[x^(n - 1) Sinh[a x], x];

(* Polynomial x Sinh^m / Cosh^m (hyperbolic analogues of circular Formulas
   381-386).  Derived by the power-reduction Sinh^2 = (Cosh 2ax - 1)/2,
   Cosh^2 = (Cosh 2ax + 1)/2, Sinh^3 = (Sinh 3ax - 3 Sinh ax)/4,
   Cosh^3 = (Cosh 3ax + 3 Cosh ax)/4; the leading x^(k+1) term flips sign
   relative to the circular case for the Sinh members (the "-1/2" constant). *)
IntegrateTable[x_ Sinh[a_. x_]^2, x_] /; FreeQ[a, x] :=
  -x^2/4 + x Sinh[2 a x]/(4 a) - Cosh[2 a x]/(8 a^2);
IntegrateTable[x_^2 Sinh[a_. x_]^2, x_] /; FreeQ[a, x] :=
  -x^3/6 + x^2 Sinh[2 a x]/(4 a) - x Cosh[2 a x]/(4 a^2) + Sinh[2 a x]/(8 a^3);
IntegrateTable[x_ Sinh[a_. x_]^3, x_] /; FreeQ[a, x] :=
  x Cosh[3 a x]/(12 a) - Sinh[3 a x]/(36 a^2) - 3 x Cosh[a x]/(4 a) + 3 Sinh[a x]/(4 a^2);
IntegrateTable[x_ Cosh[a_. x_]^2, x_] /; FreeQ[a, x] :=
  x^2/4 + x Sinh[2 a x]/(4 a) - Cosh[2 a x]/(8 a^2);
IntegrateTable[x_^2 Cosh[a_. x_]^2, x_] /; FreeQ[a, x] :=
  x^3/6 + x^2 Sinh[2 a x]/(4 a) - x Cosh[2 a x]/(4 a^2) + Sinh[2 a x]/(8 a^3);
IntegrateTable[x_ Cosh[a_. x_]^3, x_] /; FreeQ[a, x] :=
  x Sinh[3 a x]/(12 a) - Cosh[3 a x]/(36 a^2) + 3 x Sinh[a x]/(4 a) - 3 Cosh[a x]/(4 a^2);

(* Polynomial over hyperbolic powers: x Csch^n = x/Sinh^n, x Sech^n = x/Cosh^n
   (hyperbolic analogues of circular Formulas 410-413, keyed on the canonical
   Csch/Sech heads since Mathilda folds 1/Sinh^n -> Csch^n).  The recursion
   drops n by two and closes on the x Csch^2 / x Sech^2 bases (even n); odd n
   bottoms out on the non-elementary x Csch / x Sech (a dilogarithm), left
   unresolved exactly as the circular x/Sin, x/Cos are.  Sign of the recursive
   term: - for Csch (Csch^2 = Coth^2 - 1), + for Sech (Sech^2 = 1 - Tanh^2). *)
IntegrateTable[x_ Csch[a_. x_]^2, x_] /; FreeQ[a, x] :=
  -x Coth[a x]/a + Log[Sinh[a x]]/a^2;
IntegrateTable[x_ Csch[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 2 :=
  -x Coth[a x] Csch[a x]^(n - 2)/(a (n - 1)) - Csch[a x]^(n - 2)/(a^2 (n - 1) (n - 2)) - (n - 2)/(n - 1) IntegrateTable[x Csch[a x]^(n - 2), x];
IntegrateTable[x_ Sech[a_. x_]^2, x_] /; FreeQ[a, x] :=
  x Tanh[a x]/a - Log[Cosh[a x]]/a^2;
IntegrateTable[x_ Sech[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 2 :=
  x Tanh[a x] Sech[a x]^(n - 2)/(a (n - 1)) + Sech[a x]^(n - 2)/(a^2 (n - 1) (n - 2)) + (n - 2)/(n - 1) IntegrateTable[x Sech[a x]^(n - 2), x];

(* Formula 544 & 545 (generalised to argument a x) *)
IntegrateTable[Sech[a_. x_] Tanh[a_. x_], x_] /; FreeQ[a, x] := -Sech[a x]/a;
IntegrateTable[Csch[a_. x_] Coth[a_. x_], x_] /; FreeQ[a, x] := -Csch[a x]/a;

(* Formula 546 & 552 & 549 & 551 & 553 & 555 (generalised to argument a x) *)
IntegrateTable[Sinh[a_. x_]^2, x_] /; FreeQ[a, x] := Sinh[2 a x]/(4 a) - x/2;
IntegrateTable[Cosh[a_. x_]^2, x_] /; FreeQ[a, x] := Sinh[2 a x]/(4 a) + x/2;
IntegrateTable[Tanh[a_. x_]^2, x_] /; FreeQ[a, x] := x - Tanh[a x]/a;
IntegrateTable[Sech[a_. x_]^2, x_] /; FreeQ[a, x] := Tanh[a x]/a;
IntegrateTable[Coth[a_. x_]^2, x_] /; FreeQ[a, x] := x - Coth[a x]/a;
IntegrateTable[Csch[a_. x_]^2, x_] /; FreeQ[a, x] := -Coth[a x]/a;

(* Hyperbolic product base cases (mirror circular Formulas 301/304/305).
   The Formula 547 reduction below strips the Cosh power two at a time and
   bottoms out on Cosh^1 (odd n) or Cosh^0 = Sinh^m (even n); the single-Cosh
   and single-Sinh bases close the odd-n tail.  No sign flip relative to the
   circular forms: d/dx Cosh = +Sinh and d/dx Sinh = +Cosh both integrate the
   power up, whereas the circular Sin/Cos pair carries the d/dx Cos = -Sin sign
   only in the single-Sin form (Formula 304), which is why that one is negated
   there and not here. *)
IntegrateTable[Sinh[a_. x_] Cosh[a_. x_], x_] /; FreeQ[a, x] := Sinh[a x]^2/(2 a);
IntegrateTable[Sinh[a_. x_] Cosh[a_. x_]^m_, x_] /; FreeQ[{a, m}, x] && m =!= -1 :=
  Cosh[a x]^(m + 1)/((m + 1) a);
IntegrateTable[Sinh[a_. x_]^m_ Cosh[a_. x_], x_] /; FreeQ[{a, m}, x] && m =!= -1 :=
  Sinh[a x]^(m + 1)/((m + 1) a);

(* Formula 547 (generalised to argument a x) *)
IntegrateTable[Sinh[a_. x_]^m_Integer Cosh[a_. x_]^n_Integer, x_] /; FreeQ[{a, m, n}, x] && m + n =!= 0 && n > 1 :=
  (Sinh[a x]^(m + 1) Cosh[a x]^(n - 1))/((m + n) a) + (n - 1)/(m + n) IntegrateTable[Sinh[a x]^m Cosh[a x]^(n - 2), x];

(* Hyperbolic reciprocal-product reductions, keyed on the Csch / Sech HEADS.
   Mathilda canonicalises 1/Sinh -> Csch and 1/Cosh -> Sech, so the CRC
   `1/(Sinh^m Cosh^n)` form (old Formula 548) never actually matches an
   evaluated integrand — the reductions below replace it in the real head
   basis.  Signs follow 1 = Cosh^2 - Sinh^2 (the hyperbolic reduction carries a
   minus on its recursive term, unlike the circular Csc/Sec analogue's plus).
   Terminate at Csch Sech (odd m and n) or the single-head Sech^n / Csch^n. *)
IntegrateTable[Csch[a_. x_] Sech[a_. x_], x_] /; FreeQ[a, x] := Log[Tanh[a x]]/a;
IntegrateTable[Csch[a_. x_] Sech[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 1 :=
  Sech[a x]^(n - 1)/(a (n - 1)) + IntegrateTable[Csch[a x] Sech[a x]^(n - 2), x];
IntegrateTable[Csch[a_. x_]^m_Integer Sech[a_. x_]^n_Integer, x_] /; FreeQ[{a, m, n}, x] && m > 1 :=
  -(Csch[a x]^(m - 1) Sech[a x]^(n - 1))/(a (m - 1)) - (m + n - 2)/(m - 1) IntegrateTable[Csch[a x]^(m - 2) Sech[a x]^n, x];

(* Formula 550 *)
IntegrateTable[Tanh[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 1 := 
  -Tanh[a x]^(n - 1)/(a (n - 1)) + IntegrateTable[Tanh[a x]^(n - 2), x];

(* Formula 554 *)
IntegrateTable[Coth[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 1 := 
  -Coth[a x]^(n - 1)/(a (n - 1)) + IntegrateTable[Coth[a x]^(n - 2), x];

(* Sinh Reduction (missing from CRC tables) *)
IntegrateTable[Sinh[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 1 := 
  (Sinh[a x]^(n - 1) Cosh[a x])/(a n) - (n - 1)/n IntegrateTable[Sinh[a x]^(n - 2), x];

(* Cosh Reduction (missing from CRC tables) *)
IntegrateTable[Cosh[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 1 := 
  (Cosh[a x]^(n - 1) Sinh[a x])/(a n) + (n - 1)/n IntegrateTable[Cosh[a x]^(n - 2), x];

(* Sech Reduction (missing from CRC tables) *)
IntegrateTable[Sech[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 1 := 
  (Sech[a x]^(n - 2) Tanh[a x])/(a (n - 1)) + (n - 2)/(n - 1) IntegrateTable[Sech[a x]^(n - 2), x];

(* Csch Reduction (missing from CRC tables) *)
IntegrateTable[Csch[a_. x_]^n_Integer, x_] /; FreeQ[{a, n}, x] && n > 1 := 
  -(Csch[a x]^(n - 2) Coth[a x])/(a (n - 1)) - (n - 2)/(n - 1) IntegrateTable[Csch[a x]^(n - 2), x];

(* Mixed Tanh/Sech and Coth/Csch powers — hyperbolic analogues of the Tan/Sec
   and Cot/Csc reductions above; the canonical form of Sinh^m/Cosh^n and
   Cosh^m/Sinh^n.  Signs follow Sech^2 = 1 - Tanh^2 and Csch^2 = Coth^2 - 1, so
   the recursive term is + here where the circular Tan/Sec analogue is -.
   The n = 1 base (Tanh Sech, Coth Csch) also matches the Sech Tanh / Csch Coth
   rules (Formulas 544/545) with the same value. *)
IntegrateTable[Tanh[a_. x_] Sech[a_. x_]^n_., x_] /; FreeQ[{a, n}, x] && IntegerQ[n] :=
  -Sech[a x]^n/(a n);
IntegrateTable[Coth[a_. x_] Csch[a_. x_]^n_., x_] /; FreeQ[{a, n}, x] && IntegerQ[n] :=
  -Csch[a x]^n/(a n);
IntegrateTable[Tanh[a_. x_]^m_Integer Sech[a_. x_]^n_., x_] /; FreeQ[{a, m, n}, x] && IntegerQ[n] && m > 1 :=
  -Tanh[a x]^(m - 1) Sech[a x]^n/(a (m + n - 1)) + (m - 1)/(m + n - 1) IntegrateTable[Tanh[a x]^(m - 2) Sech[a x]^n, x];
IntegrateTable[Coth[a_. x_]^m_Integer Csch[a_. x_]^n_., x_] /; FreeQ[{a, m, n}, x] && IntegerQ[n] && m > 1 :=
  -Coth[a x]^(m - 1) Csch[a x]^n/(a (m + n - 1)) + (m - 1)/(m + n - 1) IntegrateTable[Coth[a x]^(m - 2) Csch[a x]^n, x];

(* Formula 556 - 558 *)
IntegrateTable[Sinh[m_. x_] Sinh[n_. x_], x_] /; FreeQ[{m, n}, x] && m^2 =!= n^2 :=
  Sinh[(m + n) x]/(2 (m + n)) - Sinh[(m - n) x]/(2 (m - n));
IntegrateTable[Cosh[m_. x_] Cosh[n_. x_], x_] /; FreeQ[{m, n}, x] && m^2 =!= n^2 := 
  Sinh[(m + n) x]/(2 (m + n)) + Sinh[(m - n) x]/(2 (m - n));
IntegrateTable[Sinh[m_. x_] Cosh[n_. x_], x_] /; FreeQ[{m, n}, x] && m^2 =!= n^2 := 
  Cosh[(m + n) x]/(2 (m + n)) + Cosh[(m - n) x]/(2 (m - n));

(* Formula 577 *)
IntegrateTable[x_^(p_ + 1) BesselJ[p_, x_], x_] /; FreeQ[p, x] := x^(p + 1) BesselJ[p + 1, x];
IntegrateTable[x_^(p_ + 1) BesselY[p_, x_], x_] /; FreeQ[p, x] := x^(p + 1) BesselY[p + 1, x];

(* Formula 578 *)
IntegrateTable[x_^(-p_ + 1) BesselJ[p_, x_], x_] /; FreeQ[p, x] := -x^(-p + 1) BesselJ[p - 1, x];
IntegrateTable[x_^(-p_ + 1) BesselY[p_, x_], x_] /; FreeQ[p, x] := -x^(-p + 1) BesselY[p - 1, x];

(* Formula 580 & 581 *)
IntegrateTable[BesselJ[1, x_], x_] := -BesselJ[0, x];
IntegrateTable[x_ BesselJ[0, x_], x_] := x BesselJ[1, x];

(* Modified Bessel I (DLMF 10.29.4): d/dx[x^p I_p] = x^p I_{p-1},
   d/dx[x^-p I_p] = x^-p I_{p+1}.  Note the second carries a '+' (vs BesselJ). *)
IntegrateTable[x_^(p_ + 1) BesselI[p_, x_], x_] /; FreeQ[p, x] := x^(p + 1) BesselI[p + 1, x];
IntegrateTable[x_^(-p_ + 1) BesselI[p_, x_], x_] /; FreeQ[p, x] := x^(-p + 1) BesselI[p - 1, x];
IntegrateTable[BesselI[1, x_], x_] := BesselI[0, x];
IntegrateTable[x_ BesselI[0, x_], x_] := x BesselI[1, x];

SetAttributes[IntegrateTable, {Protected, ReadProtected}];

(* Public wrapper.  The C dispatcher treats either head
   (Integrate`CRCTable or IntegrateTable) as a "no rule matched" signal. *)
Integrate`CRCTable[f_, x_] := IntegrateTable[f, x];
SetAttributes[Integrate`CRCTable, {Protected, ReadProtected}];

