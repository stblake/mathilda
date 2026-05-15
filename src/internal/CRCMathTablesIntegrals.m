
(* Source: CRC Standard Mathematical Tables, 31st Edition *)
(*
   Loaded lazily on first call to Integrate's CRCTable stage (see
   src/integrate.c's try_crctable / crc_lazy_load).  The table is
   addressable by users as Integrate`CRCTable[f, x]; the internal rules
   are stored on the short name IntegrateTable for readability — the
   public wrapper at the bottom of this file forwards to it.

   picocas's BeginPackage/Begin parsing inside Get is incomplete (it
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

(* Formula 43 *)
IntegrateTable[1/(c_^2 + x_^2), x_] /; FreeQ[c, x] := 1/c ArcTan[x/c];

(* Formula 44 *)
IntegrateTable[1/(c_^2 - x_^2), x_] /; FreeQ[c, x] := 1/(2 c) Log[(c + x)/(c - x)];

(* Formula 45 *)
IntegrateTable[1/(x_^2 - c_^2), x_] /; FreeQ[c, x] := 1/(2 c) Log[(x - c)/(x + c)];

(* Formula 46: Split into + and - variations *)
IntegrateTable[x_/(c_^2 + x_^2), x_] /; FreeQ[c, x] := 1/2 Log[c^2 + x^2];
IntegrateTable[x_/(c_^2 - x_^2), x_] /; FreeQ[c, x] := -1/2 Log[c^2 - x^2];

(* Formula 47: Adapted pattern to match m instead of n+1 for better matching *)
IntegrateTable[x_/(c_^2 + x_^2)^m_, x_] /; FreeQ[{c, m}, x] && m =!= 1 := -1/(2 (m - 1) (c^2 + x^2)^(m - 1));
IntegrateTable[x_/(c_^2 - x_^2)^m_, x_] /; FreeQ[{c, m}, x] && m =!= 1 := 1/(2 (m - 1) (c^2 - x^2)^(m - 1));

(* Formula 48: Split into + and - variations *)
IntegrateTable[1/(c_^2 + x_^2)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(2 c^2 (n - 1)) (x/(c^2 + x^2)^(n - 1) + (2 n - 3) IntegrateTable[1/(c^2 + x^2)^(n - 1), x]);
IntegrateTable[1/(c_^2 - x_^2)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(2 c^2 (n - 1)) (x/(c^2 - x^2)^(n - 1) + (2 n - 3) IntegrateTable[1/(c^2 - x^2)^(n - 1), x]);

(* Formula 49 *)
IntegrateTable[1/(x_^2 - c_^2)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(2 c^2 (n - 1)) (-x/(x^2 - c^2)^(n - 1) - (2 n - 3) IntegrateTable[1/(x^2 - c^2)^(n - 1), x]);

(* Formula 50 *)
IntegrateTable[x_/(x_^2 - c_^2), x_] /; FreeQ[c, x] := 1/2 Log[x^2 - c^2];

(* Formula 51: Adapted pattern to match m instead of n+1 *)
IntegrateTable[x_/(x_^2 - c_^2)^m_, x_] /; FreeQ[{c, m}, x] && m =!= 1 := -1/(2 (m - 1) (x^2 - c^2)^(m - 1));

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

(* Formula 62 *)
IntegrateTable[(a_^2 + b_^2 x_^2)^-1, x_] /; FreeQ[{a, b}, x] := 1/(a b) ArcTan[(b x)/a];

(* Formula 63 *)
IntegrateTable[x_ (a_. + b_. x_^2)^-1, x_] /; FreeQ[{a, b}, x] := 1/(2 b) Log[a + b x^2];

(* Formula 64 *)
IntegrateTable[x_^2 (a_. + b_. x_^2)^-1, x_] /; FreeQ[{a, b}, x] := x/b - a/b IntegrateTable[(a + b x^2)^-1, x];

(* Formula 65 *)
IntegrateTable[1/(a_ + b_. x_^2)^2, x_] /; FreeQ[{a, b}, x] := x/(2 a (a + b x^2)) + 1/(2 a) IntegrateTable[1/(a + b x^2), x];

(* Formula 66 *)
IntegrateTable[1/(a_^2 - b_^2 x_^2), x_] /; FreeQ[{a, b}, x] := 1/(2 a b) Log[(a + b x)/(a - b x)];

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

(* Formula 88 *)
IntegrateTable[1/(c_^3 + x_^3), x_] /; FreeQ[c, x] := 1/(6 c^2) Log[(c + x)^3/(c^3 + x^3)] + 1/(c^2 Sqrt[3]) ArcTan[(2 x - c)/(c Sqrt[3])];

IntegrateTable[1/(c_^3 - x_^3), x_] /; FreeQ[c, x] := -1/(6 c^2) Log[(c - x)^3/(c^3 - x^3)] + 1/(c^2 Sqrt[3]) ArcTan[(2 x + c)/(c Sqrt[3])];

(* Formula 89 *)
IntegrateTable[1/(c_^3 + x_^3)^2, x_] /; FreeQ[c, x] := x/(3 c^3 (c^3 + x^3)) + 2/(3 c^3) IntegrateTable[1/(c^3 + x^3), x];

IntegrateTable[1/(c_^3 - x_^3)^2, x_] /; FreeQ[c, x] := x/(3 c^3 (c^3 - x^3)) + 2/(3 c^3) IntegrateTable[1/(c^3 - x^3), x];

(* Formula 90: Mapped n+1 to n *)
IntegrateTable[1/(c_^3 + x_^3)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(3 (n - 1) c^3) (x/(c^3 + x^3)^(n - 1) + (3 n - 4) IntegrateTable[1/(c^3 + x^3)^(n - 1), x]);

IntegrateTable[1/(c_^3 - x_^3)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(3 (n - 1) c^3) (x/(c^3 - x^3)^(n - 1) + (3 n - 4) IntegrateTable[1/(c^3 - x^3)^(n - 1), x]);

(* Formula 91 *)
IntegrateTable[x_/(c_^3 + x_^3), x_] /; FreeQ[c, x] := 1/(6 c) Log[(c^3 + x^3)/(c + x)^3] + 1/(c Sqrt[3]) ArcTan[(2 x - c)/(c Sqrt[3])];

IntegrateTable[x_/(c_^3 - x_^3), x_] /; FreeQ[c, x] := 1/(6 c) Log[(c^3 - x^3)/(c - x)^3] - 1/(c Sqrt[3]) ArcTan[(2 x + c)/(c Sqrt[3])];

(* Formula 92 *)
IntegrateTable[x_/(c_^3 + x_^3)^2, x_] /; FreeQ[c, x] := x^2/(3 c^3 (c^3 + x^3)) + 1/(3 c^3) IntegrateTable[x/(c^3 + x^3), x];

IntegrateTable[x_/(c_^3 - x_^3)^2, x_] /; FreeQ[c, x] := x^2/(3 c^3 (c^3 - x^3)) + 1/(3 c^3) IntegrateTable[x/(c^3 - x^3), x];

(* Formula 93: Mapped n+1 to n *)
IntegrateTable[x_/(c_^3 + x_^3)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(3 (n - 1) c^3) (x^2/(c^3 + x^3)^(n - 1) + (3 n - 5) IntegrateTable[x/(c^3 + x^3)^(n - 1), x]);

IntegrateTable[x_/(c_^3 - x_^3)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(3 (n - 1) c^3) (x^2/(c^3 - x^3)^(n - 1) + (3 n - 5) IntegrateTable[x/(c^3 - x^3)^(n - 1), x]);

(* Formula 94 *)
IntegrateTable[x_^2/(c_^3 + x_^3), x_] /; FreeQ[c, x] := 1/3 Log[c^3 + x^3];

IntegrateTable[x_^2/(c_^3 - x_^3), x_] /; FreeQ[c, x] := -1/3 Log[c^3 - x^3];

(* Formula 95: Mapped n+1 to n *)
IntegrateTable[x_^2/(c_^3 + x_^3)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 := -1/(3 (n - 1) (c^3 + x^3)^(n - 1));

IntegrateTable[x_^2/(c_^3 - x_^3)^n_, x_] /; FreeQ[{c, n}, x] && n =!= 1 := 1/(3 (n - 1) (c^3 - x^3)^(n - 1));

(* Formula 96 *)
IntegrateTable[1/(x_ (c_^3 + x_^3)), x_] /; FreeQ[c, x] := 1/(3 c^3) Log[x^3/(c^3 + x^3)];

IntegrateTable[1/(x_ (c_^3 - x_^3)), x_] /; FreeQ[c, x] := 1/(3 c^3) Log[x^3/(c^3 - x^3)];

(* Formula 97 *)
IntegrateTable[1/(x_ (c_^3 + x_^3)^2), x_] /; FreeQ[c, x] := 1/(3 c^3 (c^3 + x^3)) + 1/(3 c^6) Log[x^3/(c^3 + x^3)];

IntegrateTable[1/(x_ (c_^3 - x_^3)^2), x_] /; FreeQ[c, x] := 1/(3 c^3 (c^3 - x^3)) + 1/(3 c^6) Log[x^3/(c^3 - x^3)];

(* Formula 98: Mapped n+1 to n *)
IntegrateTable[1/(x_ (c_^3 + x_^3)^n_), x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(3 (n - 1) c^3 (c^3 + x^3)^(n - 1)) + 1/c^3 IntegrateTable[1/(x (c^3 + x^3)^(n - 1)), x];

IntegrateTable[1/(x_ (c_^3 - x_^3)^n_), x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/(3 (n - 1) c^3 (c^3 - x^3)^(n - 1)) + 1/c^3 IntegrateTable[1/(x (c^3 - x^3)^(n - 1)), x];

(* Formula 99 *)
IntegrateTable[1/(x_^2 (c_^3 + x_^3)), x_] /; FreeQ[c, x] := -1/(c^3 x) - 1/c^3 IntegrateTable[x/(c^3 + x^3), x];

IntegrateTable[1/(x_^2 (c_^3 - x_^3)), x_] /; FreeQ[c, x] := -1/(c^3 x) + 1/c^3 IntegrateTable[x/(c^3 - x^3), x];

(* Formula 100: Mapped n+1 to n *)
IntegrateTable[1/(x_^2 (c_^3 + x_^3)^n_), x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/c^3 IntegrateTable[1/(x^2 (c^3 + x^3)^(n - 1)), x] - 1/c^3 IntegrateTable[x/(c^3 + x^3)^n, x];

IntegrateTable[1/(x_^2 (c_^3 - x_^3)^n_), x_] /; FreeQ[{c, n}, x] && n =!= 1 && IntegerQ[n] && n > 1 := 1/c^3 IntegrateTable[1/(x^2 (c^3 - x^3)^(n - 1)), x] + 1/c^3 IntegrateTable[x/(c^3 - x^3)^n, x];

(* Formula 101 *)
IntegrateTable[1/(c_^4 + x_^4), x_] /; FreeQ[c, x] := 1/(2 c^3 Sqrt[2]) (1/2 Log[(x^2 + c x Sqrt[2] + c^2)/(x^2 - c x Sqrt[2] + c^2)] + ArcTan[(c x Sqrt[2])/(c^2 - x^2)]);

(* Formula 102 *)
IntegrateTable[1/(c_^4 - x_^4), x_] /; FreeQ[c, x] := 1/(2 c^3) (1/2 Log[(c + x)/(c - x)] + ArcTan[x/c]);

(* Formula 103 *)
IntegrateTable[x_/(c_^4 + x_^4), x_] /; FreeQ[c, x] := 1/(2 c^2) ArcTan[x^2/c^2];

(* Formula 104 *)
IntegrateTable[x_/(c_^4 - x_^4), x_] /; FreeQ[c, x] := 1/(4 c^2) Log[(c^2 + x^2)/(c^2 - x^2)];

(* Formula 105 *)
IntegrateTable[x_^2/(c_^4 + x_^4), x_] /; FreeQ[c, x] := 1/(2 c Sqrt[2]) (1/2 Log[(x^2 - c x Sqrt[2] + c^2)/(x^2 + c x Sqrt[2] + c^2)] + ArcTan[(c x Sqrt[2])/(c^2 - x^2)]);

(* Formula 106 *)
IntegrateTable[x_^2/(c_^4 - x_^4), x_] /; FreeQ[c, x] := 1/(2 c) (1/2 Log[(c + x)/(c - x)] - ArcTan[x/c]);

(* Formula 107 *)
IntegrateTable[x_^3/(c_^4 + x_^4), x_] /; FreeQ[c, x] := 1/4 Log[c^4 + x^4];

IntegrateTable[x_^3/(c_^4 - x_^4), x_] /; FreeQ[c, x] := -1/4 Log[c^4 - x^4];

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

(* Formula 125 *)
IntegrateTable[x_^m_ Sqrt[a_ + b_. x_], x_] /; FreeQ[{a, b, m}, x] := 2/(b (2 m + 3)) (x^m (a + b x)^(3/2) - m a IntegrateTable[x^(m - 1) Sqrt[a + b x], x]);

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

(* Formula 131 *)
IntegrateTable[x_^m_/Sqrt[a_ + b_. x_], x_] /; FreeQ[{a, b, m}, x] := 2/((2 m + 1) b) (x^m Sqrt[a + b x] - m a IntegrateTable[x^(m - 1)/Sqrt[a + b x], x]);

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
IntegrateTable[1/((a_ + b_. x_) Sqrt[(a_ + b_. x_) (c_ + d_. x_)]), x_] /; FreeQ[{a, b, c, d}, x] && (a d - b c) d > 0 := With[{u = a + b x, v = c + d x, k = a d - b c}, 
    2/(k Sqrt[k d]) ArcTanh[Sqrt[k d u v]/(d u)]];

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

(* Formula 149 *)
IntegrateTable[Sqrt[x_^2 + a_^2], x_] /; FreeQ[a, x] := 1/2 (x Sqrt[x^2 + a^2] + a^2 Log[x + Sqrt[x^2 + a^2]]);

IntegrateTable[Sqrt[x_^2 - a_^2], x_] /; FreeQ[a, x] := 1/2 (x Sqrt[x^2 - a^2] - a^2 Log[x + Sqrt[x^2 - a^2]]);

(* Formula 150 *)
IntegrateTable[1/Sqrt[x_^2 + a_^2], x_] /; FreeQ[a, x] := Log[x + Sqrt[x^2 + a^2]];

IntegrateTable[1/Sqrt[x_^2 - a_^2], x_] /; FreeQ[a, x] := Log[x + Sqrt[x^2 - a^2]];

(* Formula 151 *)
IntegrateTable[1/(x_ Sqrt[x_^2 - a_^2]), x_] /; FreeQ[a, x] := 1/Abs[a] ArcSec[x/a];

(* Formula 152 *)
IntegrateTable[1/(x_ Sqrt[x_^2 + a_^2]), x_] /; FreeQ[a, x] := -1/a Log[(a + Sqrt[x^2 + a^2])/x];

(* Formula 153 *)
IntegrateTable[Sqrt[x_^2 + a_^2]/x_, x_] /; FreeQ[a, x] := Sqrt[x^2 + a^2] - a Log[(a + Sqrt[x^2 + a^2])/x];

(* Formula 154 *)
IntegrateTable[Sqrt[x_^2 - a_^2]/x_, x_] /; FreeQ[a, x] := Sqrt[x^2 - a^2] - Abs[a] ArcSec[x/a];

(* Formula 155 *)
IntegrateTable[x_/Sqrt[x_^2 + a_^2], x_] /; FreeQ[a, x] := Sqrt[x^2 + a^2];
IntegrateTable[x_/Sqrt[x_^2 - a_^2], x_] /; FreeQ[a, x] := Sqrt[x^2 - a^2];

(* Formula 156 *)
IntegrateTable[x_ Sqrt[x_^2 + a_^2], x_] /; FreeQ[a, x] := 1/3 (x^2 + a^2)^(3/2);
IntegrateTable[x_ Sqrt[x_^2 - a_^2], x_] /; FreeQ[a, x] := 1/3 (x^2 - a^2)^(3/2);

(* Formula 157 *)
IntegrateTable[(x_^2 + a_^2)^(3/2), x_] /; FreeQ[a, x] := 1/4 (x (x^2 + a^2)^(3/2) + (3 a^2 x)/2 Sqrt[x^2 + a^2] + (3 a^4)/2 Log[x + Sqrt[x^2 + a^2]]);

IntegrateTable[(x_^2 - a_^2)^(3/2), x_] /; FreeQ[a, x] := 1/4 (x (x^2 - a^2)^(3/2) - (3 a^2 x)/2 Sqrt[x^2 - a^2] + (3 a^4)/2 Log[x + Sqrt[x^2 - a^2]]);

(* Formula 158 *)
IntegrateTable[1/(x_^2 + a_^2)^(3/2), x_] /; FreeQ[a, x] := x/(a^2 Sqrt[x^2 + a^2]);
IntegrateTable[1/(x_^2 - a_^2)^(3/2), x_] /; FreeQ[a, x] := -x/(a^2 Sqrt[x^2 - a^2]);

(* Formula 159 *)
IntegrateTable[x_/(x_^2 + a_^2)^(3/2), x_] /; FreeQ[a, x] := -1/Sqrt[x^2 + a^2];
IntegrateTable[x_/(x_^2 - a_^2)^(3/2), x_] /; FreeQ[a, x] := -1/Sqrt[x^2 - a^2];

(* Formula 160 *)
IntegrateTable[x_ (x_^2 + a_^2)^(3/2), x_] /; FreeQ[a, x] := 1/5 (x^2 + a^2)^(5/2);
IntegrateTable[x_ (x_^2 - a_^2)^(3/2), x_] /; FreeQ[a, x] := 1/5 (x^2 - a^2)^(5/2);

(* Formula 161 *)
IntegrateTable[x_^2 Sqrt[x_^2 + a_^2], x_] /; FreeQ[a, x] := x/4 (x^2 + a^2)^(3/2) - (a^2 x)/8 Sqrt[x^2 + a^2] - a^4/8 Log[x + Sqrt[x^2 + a^2]];

IntegrateTable[x_^2 Sqrt[x_^2 - a_^2], x_] /; FreeQ[a, x] := x/4 (x^2 - a^2)^(3/2) + (a^2 x)/8 Sqrt[x^2 - a^2] - a^4/8 Log[x + Sqrt[x^2 - a^2]];

(* Formula 162 *)
IntegrateTable[x_^3 Sqrt[x_^2 + a_^2], x_] /; FreeQ[a, x] := 1/15 (3 x^2 - 2 a^2) (x^2 + a^2)^(3/2);

(* Formula 163 *)
IntegrateTable[x_^3 Sqrt[x_^2 - a_^2], x_] /; FreeQ[a, x] := 1/5 (x^2 - a^2)^(5/2) + a^2/3 (x^2 - a^2)^(3/2);

(* Formula 164 *)
IntegrateTable[x_^2/Sqrt[x_^2 + a_^2], x_] /; FreeQ[a, x] := x/2 Sqrt[x^2 + a^2] - a^2/2 Log[x + Sqrt[x^2 + a^2]];

IntegrateTable[x_^2/Sqrt[x_^2 - a_^2], x_] /; FreeQ[a, x] := x/2 Sqrt[x^2 - a^2] + a^2/2 Log[x + Sqrt[x^2 - a^2]];

(* Formula 165 *)
IntegrateTable[x_^3/Sqrt[x_^2 + a_^2], x_] /; FreeQ[a, x] := 1/3 (x^2 + a^2)^(3/2) - a^2 Sqrt[x^2 + a^2];
IntegrateTable[x_^3/Sqrt[x_^2 - a_^2], x_] /; FreeQ[a, x] := 1/3 (x^2 - a^2)^(3/2) + a^2 Sqrt[x^2 - a^2];

(* Formula 166 *)
IntegrateTable[1/(x_^2 Sqrt[x_^2 + a_^2]), x_] /; FreeQ[a, x] := -Sqrt[x^2 + a^2]/(a^2 x);
IntegrateTable[1/(x_^2 Sqrt[x_^2 - a_^2]), x_] /; FreeQ[a, x] := Sqrt[x^2 - a^2]/(a^2 x);

(* Formula 167 *)
IntegrateTable[1/(x_^3 Sqrt[x_^2 + a_^2]), x_] /; FreeQ[a, x] := -Sqrt[x^2 + a^2]/(2 a^2 x^2) + 1/(2 a^3) Log[(a + Sqrt[x^2 + a^2])/x];

(* Formula 168 *)
IntegrateTable[1/(x_^3 Sqrt[x_^2 - a_^2]), x_] /; FreeQ[a, x] := Sqrt[x^2 - a^2]/(2 a^2 x^2) + 1/(2 Abs[a]^3) ArcSec[x/a];

(* Formula 169 *)
IntegrateTable[x_^2 (x_^2 + a_^2)^(3/2), x_] /; FreeQ[a, x] := x/6 (x^2 + a^2)^(5/2) - (a^2 x)/24 (x^2 + a^2)^(3/2) - (a^4 x)/16 Sqrt[x^2 + a^2] - a^6/16 Log[x + Sqrt[x^2 + a^2]];

IntegrateTable[x_^2 (x_^2 - a_^2)^(3/2), x_] /; FreeQ[a, x] := x/6 (x^2 - a^2)^(5/2) + (a^2 x)/24 (x^2 - a^2)^(3/2) - (a^4 x)/16 Sqrt[x^2 - a^2] + a^6/16 Log[x + Sqrt[x^2 - a^2]];

(* Formula 170 *)
IntegrateTable[x_^3 (x_^2 + a_^2)^(3/2), x_] /; FreeQ[a, x] := 1/7 (x^2 + a^2)^(7/2) - a^2/5 (x^2 + a^2)^(5/2);
IntegrateTable[x_^3 (x_^2 - a_^2)^(3/2), x_] /; FreeQ[a, x] := 1/7 (x^2 - a^2)^(7/2) + a^2/5 (x^2 - a^2)^(5/2);

(* Formula 171 *)
IntegrateTable[Sqrt[x_^2 + a_^2]/x_^2, x_] /; FreeQ[a, x] := -Sqrt[x^2 + a^2]/x + Log[x + Sqrt[x^2 + a^2]];
IntegrateTable[Sqrt[x_^2 - a_^2]/x_^2, x_] /; FreeQ[a, x] := -Sqrt[x^2 - a^2]/x + Log[x + Sqrt[x^2 - a^2]];

(* Formula 172 *)
IntegrateTable[Sqrt[x_^2 + a_^2]/x_^3, x_] /; FreeQ[a, x] := -Sqrt[x^2 + a^2]/(2 x^2) - 1/(2 a) Log[(a + Sqrt[x^2 + a^2])/x];

(* Formula 173 *)
IntegrateTable[Sqrt[x_^2 - a_^2]/x_^3, x_] /; FreeQ[a, x] := -Sqrt[x^2 - a^2]/(2 x^2) + 1/(2 Abs[a]) ArcSec[x/a];

(* Formula 174 *)
IntegrateTable[Sqrt[x_^2 + a_^2]/x_^4, x_] /; FreeQ[a, x] := -(x^2 + a^2)^(3/2)/(3 a^2 x^3);
IntegrateTable[Sqrt[x_^2 - a_^2]/x_^4, x_] /; FreeQ[a, x] := (x^2 - a^2)^(3/2)/(3 a^2 x^3);

(* Formula 175 *)
IntegrateTable[x_^2/(x_^2 + a_^2)^(3/2), x_] /; FreeQ[a, x] := -x/Sqrt[x^2 + a^2] + Log[x + Sqrt[x^2 + a^2]];
IntegrateTable[x_^2/(x_^2 - a_^2)^(3/2), x_] /; FreeQ[a, x] := -x/Sqrt[x^2 - a^2] + Log[x + Sqrt[x^2 - a^2]];

(* Formula 176 *)
IntegrateTable[x_^3/(x_^2 + a_^2)^(3/2), x_] /; FreeQ[a, x] := Sqrt[x^2 + a^2] + a^2/Sqrt[x^2 + a^2];
IntegrateTable[x_^3/(x_^2 - a_^2)^(3/2), x_] /; FreeQ[a, x] := Sqrt[x^2 - a^2] - a^2/Sqrt[x^2 - a^2];

(* Formula 177 *)
IntegrateTable[1/(x_ (x_^2 + a_^2)^(3/2)), x_] /; FreeQ[a, x] := 1/(a^2 Sqrt[x^2 + a^2]) - 1/a^3 Log[(a + Sqrt[x^2 + a^2])/x];

(* Formula 178 *)
IntegrateTable[1/(x_ (x_^2 - a_^2)^(3/2)), x_] /; FreeQ[a, x] := -1/(a^2 Sqrt[x^2 - a^2]) - 1/Abs[a]^3 ArcSec[x/a];

(* Formula 179 *)
IntegrateTable[1/(x_^2 (x_^2 + a_^2)^(3/2)), x_] /; FreeQ[a, x] := -1/a^4 (Sqrt[x^2 + a^2]/x + x/Sqrt[x^2 + a^2]);
IntegrateTable[1/(x_^2 (x_^2 - a_^2)^(3/2)), x_] /; FreeQ[a, x] := -1/a^4 (Sqrt[x^2 - a^2]/x + x/Sqrt[x^2 - a^2]);

(* Formula 180 *)
IntegrateTable[1/(x_^3 (x_^2 + a_^2)^(3/2)), x_] /; FreeQ[a, x] := -1/(2 a^2 x^2 Sqrt[x^2 + a^2]) - 3/(2 a^4 Sqrt[x^2 + a^2]) + 3/(2 a^5) Log[(a + Sqrt[x^2 + a^2])/x];

(* Formula 181 *)
IntegrateTable[1/(x_^3 (x_^2 - a_^2)^(3/2)), x_] /; FreeQ[a, x] := 
  Sqrt[x^2 - a^2]/(2 a^2 x^2) - 3/(2 a^4 Sqrt[x^2 - a^2]) - 3/(2 Abs[a]^5) ArcSec[x/a];

(* Formula 182: Recursive reductions for positive powers *)
IntegrateTable[x_^m_/Sqrt[x_^2 + a_^2], x_] /; FreeQ[{a, m}, x] && m =!= 0 && IntegerQ[m] && m > 1 := x^(m - 1)/m Sqrt[x^2 + a^2] - ((m - 1) a^2)/m IntegrateTable[x^(m - 2)/Sqrt[x^2 + a^2], x];

IntegrateTable[x_^m_/Sqrt[x_^2 - a_^2], x_] /; FreeQ[{a, m}, x] && m =!= 0 && IntegerQ[m] && m > 1 := x^(m - 1)/m Sqrt[x^2 - a^2] + ((m - 1) a^2)/m IntegrateTable[x^(m - 2)/Sqrt[x^2 - a^2], x];

(* Formula 185: Recursive reductions for negative powers *)
IntegrateTable[1/(x_^m_ Sqrt[x_^2 + a_^2]), x_] /; FreeQ[{a, m}, x] && m =!= 1 && IntegerQ[m] && m > 2 := -Sqrt[x^2 + a^2]/((m - 1) a^2 x^(m - 1)) - (m - 2)/((m - 1) a^2) IntegrateTable[1/(x^(m - 2) Sqrt[x^2 + a^2]), x];

IntegrateTable[1/(x_^m_ Sqrt[x_^2 - a_^2]), x_] /; FreeQ[{a, m}, x] && m =!= 1 && IntegerQ[m] && m > 2 := Sqrt[x^2 - a^2]/((m - 1) a^2 x^(m - 1)) + (m - 2)/((m - 1) a^2) IntegrateTable[1/(x^(m - 2) Sqrt[x^2 - a^2]), x];

(* Formula 189 & 190 *)
IntegrateTable[1/((x_ - a_) Sqrt[x_^2 - a_^2]), x_] /; FreeQ[a, x] := -Sqrt[x^2 - a^2]/(a (x - a));
IntegrateTable[1/((x_ + a_) Sqrt[x_^2 - a_^2]), x_] /; FreeQ[a, x] := Sqrt[x^2 - a^2]/(a (x + a));

(* Formula 191 *)
IntegrateTable[Sqrt[a_^2 - x_^2], x_] /; FreeQ[a, x] := 
  1/2 (x Sqrt[a^2 - x^2] + a^2 ArcSin[x/Abs[a]]);

(* Formula 192 *)
IntegrateTable[1/Sqrt[a_^2 - x_^2], x_] /; FreeQ[a, x] := ArcSin[x/Abs[a]];

(* Formula 193 *)
IntegrateTable[1/(x_ Sqrt[a_^2 - x_^2]), x_] /; FreeQ[a, x] := -1/a Log[(a + Sqrt[a^2 - x^2])/x];

(* Formula 194 *)
IntegrateTable[Sqrt[a_^2 - x_^2]/x_, x_] /; FreeQ[a, x] := Sqrt[a^2 - x^2] - a Log[(a + Sqrt[a^2 - x^2])/x];

(* Formula 195 & 196 *)
IntegrateTable[x_/Sqrt[a_^2 - x_^2], x_] /; FreeQ[a, x] := -Sqrt[a^2 - x^2];
IntegrateTable[x_ Sqrt[a_^2 - x_^2], x_] /; FreeQ[a, x] := -1/3 (a^2 - x^2)^(3/2);

(* Formula 197 *)
IntegrateTable[(a_^2 - x_^2)^(3/2), x_] /; FreeQ[a, x] := 
  1/4 (x (a^2 - x^2)^(3/2) + (3 a^2 x)/2 Sqrt[a^2 - x^2] + (3 a^4)/2 ArcSin[x/Abs[a]]);

(* Formula 198 & 199 *)
IntegrateTable[1/(a_^2 - x_^2)^(3/2), x_] /; FreeQ[a, x] := x/(a^2 Sqrt[a^2 - x^2]);
IntegrateTable[x_/(a_^2 - x_^2)^(3/2), x_] /; FreeQ[a, x] := 1/Sqrt[a^2 - x^2];

(* Formula 200 *)
IntegrateTable[x_ (a_^2 - x_^2)^(3/2), x_] /; FreeQ[a, x] := -1/5 (a^2 - x^2)^(5/2);

(* Formula 201 *)
IntegrateTable[x_^2 Sqrt[a_^2 - x_^2], x_] /; FreeQ[a, x] := 
  -x/4 (a^2 - x^2)^(3/2) + a^2/8 (x Sqrt[a^2 - x^2] + a^2 ArcSin[x/Abs[a]]);

(* Formula 202 *)
IntegrateTable[x_^3 Sqrt[a_^2 - x_^2], x_] /; FreeQ[a, x] := (-1/5 x^2 - 2/15 a^2) (a^2 - x^2)^(3/2);

(* Formula 203 *)
IntegrateTable[x_^2 (a_^2 - x_^2)^(3/2), x_] /; FreeQ[a, x] := 
  -1/6 x (a^2 - x^2)^(5/2) + (a^2 x)/24 (a^2 - x^2)^(3/2) + (a^4 x)/16 Sqrt[a^2 - x^2] + a^6/16 ArcSin[x/Abs[a]];

(* Formula 204 *)
IntegrateTable[x_^3 (a_^2 - x_^2)^(3/2), x_] /; FreeQ[a, x] := 1/7 (a^2 - x^2)^(7/2) - a^2/5 (a^2 - x^2)^(5/2);

(* Formula 205 *)
IntegrateTable[x_^2/Sqrt[a_^2 - x_^2], x_] /; FreeQ[a, x] := 
  -x/2 Sqrt[a^2 - x^2] + a^2/2 ArcSin[x/Abs[a]];

(* Formula 206 *)
IntegrateTable[1/(x_^2 Sqrt[a_^2 - x_^2]), x_] /; FreeQ[a, x] := -Sqrt[a^2 - x^2]/(a^2 x);

(* Formula 207 *)
IntegrateTable[Sqrt[a_^2 - x_^2]/x_^2, x_] /; FreeQ[a, x] := -Sqrt[a^2 - x^2]/x - ArcSin[x/Abs[a]];

(* Formula 208 *)
IntegrateTable[Sqrt[a_^2 - x_^2]/x_^3, x_] /; FreeQ[a, x] := 
  -Sqrt[a^2 - x^2]/(2 x^2) + 1/(2 a) Log[(a + Sqrt[a^2 - x^2])/x];

(* Formula 209 *)
IntegrateTable[Sqrt[a_^2 - x_^2]/x_^4, x_] /; FreeQ[a, x] := -(a^2 - x^2)^(3/2)/(3 a^2 x^3);

(* Formula 210 *)
IntegrateTable[x_^2/(a_^2 - x_^2)^(3/2), x_] /; FreeQ[a, x] := x/Sqrt[a^2 - x^2] - ArcSin[x/Abs[a]];

(* Formula 211 *)
IntegrateTable[x_^3/Sqrt[a_^2 - x_^2], x_] /; FreeQ[a, x] := -2/3 (a^2 - x^2)^(3/2) - x^2 Sqrt[a^2 - x^2];

(* Formula 212 *)
IntegrateTable[x_^3/(a_^2 - x_^2)^(3/2), x_] /; FreeQ[a, x] := 2 Sqrt[a^2 - x^2] + x^2/Sqrt[a^2 - x^2];

(* Formula 213 *)
IntegrateTable[1/(x_^3 Sqrt[a_^2 - x_^2]), x_] /; FreeQ[a, x] := 
  -Sqrt[a^2 - x^2]/(2 a^2 x^2) - 1/(2 a^3) Log[(a + Sqrt[a^2 - x^2])/x];

(* Formula 214 *)
IntegrateTable[1/(x_ (a_^2 - x_^2)^(3/2)), x_] /; FreeQ[a, x] := 
  1/(a^2 Sqrt[a^2 - x^2]) - 1/a^3 Log[(a + Sqrt[a^2 - x^2])/x];

(* Formula 215 *)
IntegrateTable[1/(x_^2 (a_^2 - x_^2)^(3/2)), x_] /; FreeQ[a, x] := 
  1/a^4 (-Sqrt[a^2 - x^2]/x + x/Sqrt[a^2 - x^2]);

(* Formula 216 *)
IntegrateTable[1/(x_^3 (a_^2 - x_^2)^(3/2)), x_] /; FreeQ[a, x] := 
  (3 x^2 - 2 a^2)/(2 a^4 Sqrt[a^2 - x^2]) - 3/(2 a^5) Log[(a + Sqrt[a^2 - x^2])/x];

(* Formula 217 & 220: General Reduction Rules *)
IntegrateTable[x_^m_/Sqrt[a_^2 - x_^2], x_] /; FreeQ[{a, m}, x] && m =!= 0 && IntegerQ[m] && m > 1 :=
  -x^(m - 1)/m Sqrt[a^2 - x^2] + ((m - 1) a^2)/m IntegrateTable[x^(m - 2)/Sqrt[a^2 - x^2], x];
IntegrateTable[1/(x_^m_ Sqrt[a_^2 - x_^2]), x_] /; FreeQ[{a, m}, x] && m =!= 1 && IntegerQ[m] && m > 2 :=
  -Sqrt[a^2 - x^2]/((m - 1) a^2 x^(m - 1)) + (m - 2)/((m - 1) a^2) IntegrateTable[1/(x^(m - 2) Sqrt[a^2 - x^2]), x];

(* Formula 223 *)
IntegrateTable[1/((b_^2 - x_^2) Sqrt[a_^2 - x_^2]), x_] /; FreeQ[{a, b}, x] && a^2 > b^2 := 
  1/(2 b Sqrt[a^2 - b^2]) Log[(b Sqrt[a^2 - x^2] + x Sqrt[a^2 - b^2])/(b Sqrt[a^2 - x^2] - x Sqrt[a^2 - b^2])];
IntegrateTable[1/((b_^2 - x_^2) Sqrt[a_^2 - x_^2]), x_] /; FreeQ[{a, b}, x] && b^2 > a^2 := 
  1/(b Sqrt[b^2 - a^2]) ArcTan[(x Sqrt[b^2 - a^2])/(b Sqrt[a^2 - x^2])];

(* Formula 224 *)
IntegrateTable[1/((b_^2 + x_^2) Sqrt[a_^2 - x_^2]), x_] /; FreeQ[{a, b}, x] := 
  1/(b Sqrt[a^2 + b^2]) ArcTan[(x Sqrt[a^2 + b^2])/(b Sqrt[a^2 - x^2])];

(* Formula 225 *)
IntegrateTable[Sqrt[a_^2 - x_^2]/(b_^2 + x_^2), x_] /; FreeQ[{a, b}, x] := 
  Sqrt[a^2 + b^2]/Abs[b] ArcSin[(x Sqrt[a^2 + b^2])/(Abs[a] Sqrt[x^2 + b^2])] - ArcSin[x/Abs[a]];

(* Formula 226 *)
IntegrateTable[1/Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c}, x] && c > 0 := 
  1/Sqrt[c] Log[2 Sqrt[c (a + b x + c x^2)] + 2 c x + b];
IntegrateTable[1/Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c}, x] && c < 0 := 
  -1/Sqrt[-c] ArcSin[(2 c x + b)/Sqrt[b^2 - 4 a c]];

(* Formula 227 & 228 *)
IntegrateTable[1/(a_ + b_. x_ + c_. x_^2)^(3/2), x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2, q = 4 a c - b^2}, (2 (2 c x + b))/(q Sqrt[X])];
IntegrateTable[1/(a_ + b_. x_ + c_. x_^2)^(5/2), x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2, q = 4 a c - b^2, k = (4 c)/q}, 
    (2 (2 c x + b))/(3 q Sqrt[X]) (1/X + 2 k)
  ];

(* Formula 229: Generic reduction for n (mapped from texts n) *)
IntegrateTable[1/((a_ + b_. x_ + c_. x_^2)^n_ Sqrt[a_ + b_. x_ + c_. x_^2]), x_] /; FreeQ[{a, b, c, n}, x] && n =!= 0 && IntegerQ[n] && n > 0 :=
  With[{X = a + b x + c x^2, q = 4 a c - b^2, k = (4 c)/q},
    (2 (2 c x + b) Sqrt[X])/((2 n + 1) q X^(n + 1)) + (2 k n)/(2 n + 1) IntegrateTable[1/(X^(n - 1) Sqrt[X]), x]
  ];

(* Formula 230 - 232 *)
IntegrateTable[Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2, q = 4 a c - b^2, k = (4 c)/q}, 
    ((2 c x + b) Sqrt[X])/(4 c) + 1/(2 k) IntegrateTable[1/Sqrt[X], x]
  ];
IntegrateTable[(a_ + b_. x_ + c_. x_^2)^(3/2), x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2, q = 4 a c - b^2, k = (4 c)/q}, 
    ((2 c x + b) Sqrt[X])/(8 c) (X + 3/(2 k)) + 3/(8 k^2) IntegrateTable[1/Sqrt[X], x]
  ];
IntegrateTable[(a_ + b_. x_ + c_. x_^2)^(5/2), x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2, q = 4 a c - b^2, k = (4 c)/q}, 
    ((2 c x + b) Sqrt[X])/(12 c) (X^2 + (5 X)/(4 k) + 15/(8 k^2)) + 5/(16 k^3) IntegrateTable[1/Sqrt[X], x]
  ];

(* Formula 233: General reduction *)
IntegrateTable[(a_ + b_. x_ + c_. x_^2)^n_ Sqrt[a_ + b_. x_ + c_. x_^2], x_] /; FreeQ[{a, b, c, n}, x] && IntegerQ[n] && n > 0 :=
  With[{X = a + b x + c x^2, q = 4 a c - b^2, k = (4 c)/q},
    ((2 c x + b) X^n Sqrt[X])/(4 (n + 1) c) + (2 n + 1)/(2 (n + 1) k) IntegrateTable[X^(n - 1) Sqrt[X], x]
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

(* Formula 242 - 245 *)
IntegrateTable[x_ (a_ + b_. x_ + c_. x_^2)^(3/2), x_] /; FreeQ[{a, b, c}, x] := 
  With[{X = a + b x + c x^2, q = 4 a c - b^2, k = (4 c)/q}, 
    (X Sqrt[X])/(3 c) - (b (2 c x + b))/(8 c^2) Sqrt[X] - b/(4 c k) IntegrateTable[1/Sqrt[X], x]
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

(* Formula 270 *)
IntegrateTable[1/(x_ Sqrt[x_^n_ + a_^2]), x_] /; FreeQ[{a, n}, x] := 
  -2/(n a) Log[(a + Sqrt[x^n + a^2])/Sqrt[x^n]];

(* Formula 271 *)
IntegrateTable[1/(x_ Sqrt[x_^n_ - a_^2]), x_] /; FreeQ[{a, n}, x] := 
  -2/(n a) ArcSin[a/Sqrt[x^n]];

(* Formula 272 *)
IntegrateTable[Sqrt[x_/(a_^3 - x_^3)], x_] /; FreeQ[a, x] := 
  2/3 ArcSin[(x/a)^(3/2)];

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
IntegrateTable[1/(a_^2 Cos[x_]^2 + b_^2 Sin[x_]^2), x_] /; FreeQ[{a, b}, x] := 
  1/(a b) ArcTan[(b Tan[x])/a];

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
IntegrateTable[1/(a_^2 + b_^2 Sin[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] := 
  1/(a c Sqrt[a^2 + b^2]) ArcTan[(Sqrt[a^2 + b^2] Tan[c x])/a];
IntegrateTable[1/(a_^2 - b_^2 Sin[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] && a^2 > b^2 := 
  1/(a c Sqrt[a^2 - b^2]) ArcTan[(Sqrt[a^2 - b^2] Tan[c x])/a];
IntegrateTable[1/(a_^2 - b_^2 Sin[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] && b^2 > a^2 := 
  1/(2 a c Sqrt[b^2 - a^2]) Log[(Sqrt[b^2 - a^2] Tan[c x] + a)/(Sqrt[b^2 - a^2] Tan[c x] - a)];

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

(* Formula 356 *)
IntegrateTable[1/(a_^2 + b_^2 - 2 a_ b_ Cos[c_. x_]), x_] /; FreeQ[{a, b, c}, x] := 
  2/(c (a^2 - b^2)) ArcTan[((a + b) Tan[(c x)/2])/(a - b)];

(* Formula 357 & 358 *)
IntegrateTable[1/(a_^2 + b_^2 Cos[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] := 
  1/(a c Sqrt[a^2 + b^2]) ArcTan[(a Tan[c x])/Sqrt[a^2 + b^2]];
IntegrateTable[1/(a_^2 - b_^2 Cos[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] && a^2 > b^2 := 
  1/(a c Sqrt[a^2 - b^2]) ArcTan[(a Tan[c x])/Sqrt[a^2 - b^2]];
IntegrateTable[1/(a_^2 - b_^2 Cos[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] && b^2 > a^2 := 
  1/(2 a c Sqrt[b^2 - a^2]) Log[(a Tan[c x] - Sqrt[b^2 - a^2])/(a Tan[c x] + Sqrt[b^2 - a^2])];

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
IntegrateTable[1/(a_^2 Cos[c_. x_]^2 - b_^2 Sin[c_. x_]^2), x_] /; FreeQ[{a, b, c}, x] := 
  1/(2 a b c) Log[(b Tan[c x] + a)/(b Tan[c x] - a)];

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

(* Formula 396: Primary Branch *)
IntegrateTable[Sqrt[1 + Sin[x_]], x_] := 2 (Sin[x/2] - Cos[x/2]);

(* Formula 397: branch-correct.  The classical form 2 (Sin[x/2] + Cos[x/2])
   differentiates to Cos[x/2] - Sin[x/2], matching Sqrt[1 - Sin[x]] only
   when Cos[x/2] - Sin[x/2] >= 0.  Using u = Pi/2 - x and the branch-correct
   Formula 394 yields the form below, whose derivative equals
   Sqrt[1 - Sin[x]] on every branch. *)
IntegrateTable[Sqrt[1 - Sin[x_]], x_] := 2 Tan[Pi/4 + x/2] Sqrt[1 - Sin[x]];

(* Formula 398 *)
IntegrateTable[1/Sqrt[1 - Cos[x_]], x_] := Sqrt[2] Log[Tan[x/4]];

(* Formula 399 *)
IntegrateTable[1/Sqrt[1 + Cos[x_]], x_] := Sqrt[2] Log[Tan[(x + Pi)/4]];

(* Formula 400 *)
IntegrateTable[1/Sqrt[1 - Sin[x_]], x_] := Sqrt[2] Log[Tan[x/4 - Pi/8]];

(* Formula 401: Primary Branch *)
IntegrateTable[1/Sqrt[1 + Sin[x_]], x_] := 
  Sqrt[2] Log[Tan[x/4 + Pi/8]];

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

(* Formula 414 *)
IntegrateTable[Sin[a_. x_]/Sqrt[1 + b_^2 Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] := 
  -1/(a b) ArcSin[(b Cos[a x])/Sqrt[1 + b^2]];

(* Formula 415 *)
IntegrateTable[Sin[a_. x_]/Sqrt[1 - b_^2 Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] := 
  -1/(a b) Log[b Cos[a x] + Sqrt[1 - b^2 Sin[a x]^2]];

(* Formula 416 *)
IntegrateTable[Sin[a_. x_] Sqrt[1 + b_^2 Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] := 
  -(Cos[a x]/(2 a)) Sqrt[1 + b^2 Sin[a x]^2] - (1 + b^2)/(2 a b) ArcSin[(b Cos[a x])/Sqrt[1 + b^2]];

(* Formula 417 *)
IntegrateTable[Sin[a_. x_] Sqrt[1 - b_^2 Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] := 
  -(Cos[a x]/(2 a)) Sqrt[1 - b^2 Sin[a x]^2] - (1 - b^2)/(2 a b) Log[b Cos[a x] + Sqrt[1 - b^2 Sin[a x]^2]];

(* Formula 418 *)
IntegrateTable[Cos[a_. x_]/Sqrt[1 + b_^2 Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] := 
  1/(a b) Log[b Sin[a x] + Sqrt[1 + b^2 Sin[a x]^2]];

(* Formula 419 *)
IntegrateTable[Cos[a_. x_]/Sqrt[1 - b_^2 Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] := 
  1/(a b) ArcSin[b Sin[a x]];

(* Formula 420 *)
IntegrateTable[Cos[a_. x_] Sqrt[1 + b_^2 Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] := 
  (Sin[a x]/(2 a)) Sqrt[1 + b^2 Sin[a x]^2] + 1/(2 a b) Log[b Sin[a x] + Sqrt[1 + b^2 Sin[a x]^2]];

(* Formula 421 *)
IntegrateTable[Cos[a_. x_] Sqrt[1 - b_^2 Sin[a_. x_]^2], x_] /; FreeQ[{a, b}, x] := 
  (Sin[a x]/(2 a)) Sqrt[1 - b^2 Sin[a x]^2] + 1/(2 a b) ArcSin[b Sin[a x]];

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

(* Formula 453 *)
IntegrateTable[ArcSin[a_. x_]/Sqrt[1 - a_^2 x_^2], x_] /; FreeQ[a, x] := 
  ArcSin[a x]^2/(2 a);

(* Formula 454: Mapped explicitly from the integral recurrence relation *)
IntegrateTable[x_^n_ ArcSin[a_. x_]/Sqrt[1 - a_^2 x_^2], x_] /; FreeQ[{a, n}, x] && n =!= 0 && IntegerQ[n] && n > 1 :=
  -(x^(n - 1)/(n a^2)) Sqrt[1 - a^2 x^2] ArcSin[a x] + x^n/(n^2 a) + (n - 1)/(n a^2) IntegrateTable[(x^(n - 2) ArcSin[a x])/Sqrt[1 - a^2 x^2], x];

(* Formula 455 *)
IntegrateTable[ArcCos[a_. x_]/Sqrt[1 - a_^2 x_^2], x_] /; FreeQ[a, x] := 
  -ArcCos[a x]^2/(2 a);

(* Formula 456 *)
IntegrateTable[x_^n_ ArcCos[a_. x_]/Sqrt[1 - a_^2 x_^2], x_] /; FreeQ[{a, n}, x] && n =!= 0 && IntegerQ[n] && n > 1 :=
  -(x^(n - 1)/(n a^2)) Sqrt[1 - a^2 x^2] ArcCos[a x] - x^n/(n^2 a) + (n - 1)/(n a^2) IntegrateTable[(x^(n - 2) ArcCos[a x])/Sqrt[1 - a^2 x^2], x];

(* Formula 457 *)
IntegrateTable[ArcTan[a_. x_]/(1 + a_^2 x_^2), x_] /; FreeQ[a, x] := 
  ArcTan[a x]^2/(2 a);

(* Formula 458 *)
IntegrateTable[ArcCot[a_. x_]/(1 + a_^2 x_^2), x_] /; FreeQ[a, x] := 
  -ArcCot[a x]^2/(2 a);

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

(* Formula 488 *)
IntegrateTable[Log[x_^2 + a_^2], x_] /; FreeQ[a, x] := 
  x Log[x^2 + a^2] - 2 x + 2 a ArcTan[x/a];

(* Formula 489 *)
IntegrateTable[Log[x_^2 - a_^2], x_] /; FreeQ[a, x] := 
  x Log[x^2 - a^2] - 2 x + a Log[(x + a)/(x - a)];

(* Formula 490 *)
IntegrateTable[x_ Log[x_^2 + a_^2], x_] /; FreeQ[a, x] := 
  1/2 (x^2 + a^2) Log[x^2 + a^2] - x^2/2;

(* Formula 491 *)
IntegrateTable[Log[x_ + Sqrt[x_^2 + a_^2]], x_] /; FreeQ[a, x] := 
  x Log[x + Sqrt[x^2 + a^2]] - Sqrt[x^2 + a^2];
IntegrateTable[Log[x_ + Sqrt[x_^2 - a_^2]], x_] /; FreeQ[a, x] := 
  x Log[x + Sqrt[x^2 - a^2]] - Sqrt[x^2 - a^2];

(* Formula 492 *)
IntegrateTable[x_ Log[x_ + Sqrt[x_^2 + a_^2]], x_] /; FreeQ[a, x] := 
  (x^2/2 + a^2/4) Log[x + Sqrt[x^2 + a^2]] - (x Sqrt[x^2 + a^2])/4;
IntegrateTable[x_ Log[x_ + Sqrt[x_^2 - a_^2]], x_] /; FreeQ[a, x] := 
  (x^2/2 - a^2/4) Log[x + Sqrt[x^2 - a^2]] - (x Sqrt[x^2 - a^2])/4;

(* Formula 493 *)
IntegrateTable[x_^m_ Log[x_ + Sqrt[x_^2 + a_^2]], x_] /; FreeQ[{a, m}, x] := 
  x^(m + 1)/(m + 1) Log[x + Sqrt[x^2 + a^2]] - 1/(m + 1) IntegrateTable[x^(m + 1)/Sqrt[x^2 + a^2], x];
IntegrateTable[x_^m_ Log[x_ + Sqrt[x_^2 - a_^2]], x_] /; FreeQ[{a, m}, x] := 
  x^(m + 1)/(m + 1) Log[x + Sqrt[x^2 - a^2]] - 1/(m + 1) IntegrateTable[x^(m + 1)/Sqrt[x^2 - a^2], x];

(* Formula 494 *)
IntegrateTable[Log[x_ + Sqrt[x_^2 + a_^2]]/x_^2, x_] /; FreeQ[a, x] := 
  -Log[x + Sqrt[x^2 + a^2]]/x - 1/a Log[(a + Sqrt[x^2 + a^2])/x];

(* Formula 495 *)
IntegrateTable[Log[x_ + Sqrt[x_^2 - a_^2]]/x_^2, x_] /; FreeQ[a, x] := 
  -Log[x + Sqrt[x^2 - a^2]]/x + 1/Abs[a] ArcSec[x/a];

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

(* Formula 528: Form 1 *)
IntegrateTable[E^(a_. x_) Cos[x_]^m_ Sin[x_]^n_, x_] /; FreeQ[{a, m, n}, x] && m + n;

(* Formula 534 - 539 *)
IntegrateTable[Sinh[x_], x_] := Cosh[x];
IntegrateTable[Cosh[x_], x_] := Sinh[x];
IntegrateTable[Tanh[x_], x_] := Log[Cosh[x]];
IntegrateTable[Coth[x_], x_] := Log[Sinh[x]];
IntegrateTable[Sech[x_], x_] := ArcTan[Sinh[x]];
IntegrateTable[Csch[x_], x_] := Log[Tanh[x/2]];

(* Formula 540 & 542 *)
IntegrateTable[x_ Sinh[x_], x_] := x Cosh[x] - Sinh[x];
IntegrateTable[x_ Cosh[x_], x_] := x Sinh[x] - Cosh[x];

(* Formula 541 & 543 *)
IntegrateTable[x_^n_ Sinh[x_], x_] /; FreeQ[n, x] && n > 0 := 
  x^n Cosh[x] - n IntegrateTable[x^(n - 1) Cosh[x], x];
IntegrateTable[x_^n_ Cosh[x_], x_] /; FreeQ[n, x] && n > 0 := 
  x^n Sinh[x] - n IntegrateTable[x^(n - 1) Sinh[x], x];

(* Formula 544 & 545 *)
IntegrateTable[Sech[x_] Tanh[x_], x_] := -Sech[x];
IntegrateTable[Csch[x_] Coth[x_], x_] := -Csch[x];

(* Formula 546 & 552 & 549 & 551 & 553 & 555 *)
IntegrateTable[Sinh[x_]^2, x_] := Sinh[2 x]/4 - x/2;
IntegrateTable[Cosh[x_]^2, x_] := Sinh[2 x]/4 + x/2;
IntegrateTable[Tanh[x_]^2, x_] := x - Tanh[x];
IntegrateTable[Sech[x_]^2, x_] := Tanh[x];
IntegrateTable[Coth[x_]^2, x_] := x - Coth[x];
IntegrateTable[Csch[x_]^2, x_] := -Coth[x];

(* Formula 547 *)
IntegrateTable[Sinh[x_]^m_Integer Cosh[x_]^n_Integer, x_] /; FreeQ[{m, n}, x] && m + n =!= 0 && n > 1 :=
  (Sinh[x]^(m + 1) Cosh[x]^(n - 1))/(m + n) + (n - 1)/(m + n) IntegrateTable[Sinh[x]^m Cosh[x]^(n - 2), x];

(* Formula 548 *)
IntegrateTable[1/(Sinh[x_]^m_Integer Cosh[x_]^n_Integer), x_] /; FreeQ[{m, n}, x] && m =!= 1 && m > 1 :=
  -1/((m - 1) Sinh[x]^(m - 1) Cosh[x]^(n - 1)) - (m + n - 2)/(m - 1) IntegrateTable[1/(Sinh[x]^(m - 2) Cosh[x]^n), x];

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

SetAttributes[IntegrateTable, {Protected, ReadProtected}];

(* Public wrapper.  The C dispatcher treats either head
   (Integrate`CRCTable or IntegrateTable) as a "no rule matched" signal. *)
Integrate`CRCTable[f_, x_] := IntegrateTable[f, x];
SetAttributes[Integrate`CRCTable, {Protected, ReadProtected}];

