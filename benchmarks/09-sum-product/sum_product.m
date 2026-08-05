(* ==========================================================================
   Experiment 09 -- Symbolic summation and products
   ==========================================================================
   WHAT IT MEASURES.  src/sum/ (16 files: Gosper, hypergeometric, Euler-
   Maclaurin, zeta series, telescoping) and src/product/ (14 files).

   Symbolic summation is where a CAS most visibly either has an algorithm or
   does not: Gosper's algorithm either finds the hypergeometric antidifference
   or proves none exists.  There is no slow-but-correct middle, so these rows
   split cleanly into "fast" and "absent" -- exactly the distinction ABSENT.md
   exists to keep separate from speed.

   Checks are exact closed-form VALUES, substituted at a small n.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"Sum", "Product", "Zeta", "Binomial", "Factorial",
         "HarmonicNumber", "Infinity"}];

(* 1. A polynomial sum with a Faulhaber closed form. *)
bench["Sum k^5 to n, closed form", Sum[k^5, {k, 1, n}];];
check["Sum k^5 to n, closed form", Sum[k^5, {k, 1, n}] /. n -> 10];

(* 2. A convergent zeta value. *)
bench["Sum 1/k^2 to Infinity", Sum[1/k^2, {k, 1, Infinity}];];
check["Sum 1/k^2 to Infinity", Round[10^6 N[Sum[1/k^2, {k, 1, Infinity}]]]];

(* 3. A geometric sum with a symbolic ratio. *)
bench["Sum r^k to n, symbolic ratio", Sum[r^k, {k, 0, n}];];
check["Sum r^k to n, symbolic ratio", Sum[r^k, {k, 0, n}] /. {r -> 3, n -> 6}];

(* 4. The binomial row sum -- 2^n, via Gosper or a known identity. *)
bench["Sum Binomial[n,k] over k", Sum[Binomial[n, k], {k, 0, n}];];
check["Sum Binomial[n,k] over k", Sum[Binomial[n, k], {k, 0, n}] /. n -> 12];

(* 5. A telescoping rational sum. *)
bench["Sum 1/(k(k+1)) to n", Sum[1/(k (k + 1)), {k, 1, n}];];
check["Sum 1/(k(k+1)) to n",
  Round[10^6 N[Sum[1/(k (k + 1)), {k, 1, n}] /. n -> 10]]];

(* 6. A product with a Gamma closed form. *)
bench["Product (1+1/k) to n", Product[1 + 1/k, {k, 1, n}];];
check["Product (1+1/k) to n", Product[1 + 1/k, {k, 1, n}] /. n -> 10];
