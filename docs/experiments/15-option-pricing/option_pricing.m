(* ==========================================================================
   Experiment 15 -- Option pricing: trees, PDEs, and the positive part
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.

       Mathilda      ./Mathilda -file option_pricing.m
       Mathematica   wolframscript -file option_pricing.m
       Python        python3 option_pricing.py

   WHAT IT MEASURES.  Two things no other benchmark does.

   1. A WORKING VECTOR THAT SHRINKS.  A binomial tree starts with 4001 nodes
      and ends with 1, so its last few hundred iterations run BELOW the
      packing threshold.  Every packed-array benchmark so far has been
      comfortably above the threshold or comfortably below it; none has
      CROSSED it inside a loop.

   2. AN EARLY-EXERCISE PROJECTION.  An American option takes
      max(continuation, payoff) at every node of every time step.  That is the
      positive part -- the same operation as a rectified linear unit, and the
      most common nonlinearity in numerical code.

   TWO PRICERS, ON PURPOSE.  A 4000-step Cox-Ross-Rubinstein tree and a
   1000 x 25000 explicit finite-difference solver in log-price.  They are
   independent discretisations of the same problem and agree to about 0.6%,
   which is a much stronger statement than either agreeing with itself across
   three systems: it says the two programs are pricing the same option.

   THE FD STABILITY CONDITION IS LOAD-BEARING.  An explicit scheme needs
   sigma^2 dt / dx^2 <= 1/2.  With 1000 space points over [-1.2, 1.2] and
   25000 time steps that ratio is 0.433.  The first draft of this benchmark
   used 2000 x 7000, where it is 6.19; the scheme diverged, every system
   returned NaN or a plausible-looking wrong number, and only the
   cross-system value check caught it.
   ========================================================================== *)

(* ---- shared reporting helpers (identical in every experiment file) ------- *)

SetAttributes[bench, HoldRest];

bench[label_String, expr_] := Module[{ts},
  expr;
  ts = Table[First[AbsoluteTiming[expr]], {3}];
  Print[StringPadRight[label, 52], ToString[Round[1000. Min[ts], 0.001]], " ms"]
];

check[label_String, value_] :=
  Print[StringPadRight[label, 52], "check = ", value];

(* The positive part.  Ramp is the idiomatic spelling and the one this
   experiment caused to exist in Mathilda; see README.md.  Kept as a named
   helper so the alternative spellings can be swapped in and measured. *)
relu[zz_] := Ramp[zz];

(* ---- contract ----------------------------------------------------------- *)

btn   = 4000;      (* tree steps                                            *)
btt   = 1.;        (* years to expiry                                       *)
btr   = 0.03;      (* risk-free rate                                        *)
btsig = 0.25;      (* volatility                                            *)
btk   = 100.;      (* strike                                                *)
bts0  = 100.;      (* spot                                                  *)

(* ---- 1. Cox-Ross-Rubinstein binomial tree, American put ----------------- *)

btdt = btt/btn;
btu  = Exp[btsig Sqrt[btdt]];                        (* up factor           *)
btp  = (Exp[btr btdt] - 1./btu)/(btu - 1./btu);      (* risk-neutral prob   *)
btq  = 1. - btp;
btdc = Exp[-btr btdt];                               (* one-step discount   *)

(* Backward induction.  At level k the spot lattice is S0 u^(2j-k); going one
   level back multiplies the surviving nodes by u, so the spot vector is
   maintained by `btu Most[s]` rather than rebuilt.  The vector shrinks by one
   element per iteration, which is the point of this row. *)
amtree[] := Module[{v, s, k},
  s = Table[bts0 btu^(2. j - btn), {j, 0, btn}];
  v = relu[btk - s];                                 (* terminal payoff     *)
  k = btn;
  While[k > 0,
    s = btu Most[s];
    v = btdc (btp Rest[v] + btq Most[v]);            (* continuation value  *)
    v = MapThread[Max, {v, relu[btk - s]}];          (* early exercise      *)
    k = k - 1];
  First[v]];

(* ---- 2. explicit finite difference in log-price, American put ----------- *)

fdm  = 1000;       (* space points                                          *)
fdn  = 25000;      (* time steps                                            *)
fdxl = 1.2;        (* log-price half-width; S ranges over 100 e^(+-1.2)     *)

fddx = 2. fdxl/(fdm - 1);
fddt = 1./fdn;
fdx   = Table[-fdxl + (j - 1) fddx, {j, 1, fdm}];
fdpay = relu[btk - bts0 Exp[fdx]];                   (* payoff, once        *)

(* Constant coefficients: the log transform removes the S-dependence, so the
   three stencil weights are scalars and the sweep is one shifted triple. *)
fda = 0.5 btsig^2 fddt/fddx^2 - 0.5 (btr - 0.5 btsig^2) fddt/fddx;
fdc = 0.5 btsig^2 fddt/fddx^2 + 0.5 (btr - 0.5 btsig^2) fddt/fddx;
fdb = 1. - btsig^2 fddt/fddx^2 - btr fddt;

fdlo  = btk - bts0 Exp[-fdxl];                       (* deep in the money   *)
fdmid = Quotient[fdm + 1, 2];                        (* the at-the-money x  *)

(* Dirichlet boundaries by construction: the stencil is applied to the
   interior only and the two edges are re-attached with Join.  That Join is
   the reason this row exists -- see README.md. *)
amfd[] := Module[{v, k, vi},
  v = fdpay;
  k = 0;
  While[k < fdn,
    vi = fda Most[Most[v]] + fdb Take[v, {2, -2}] + fdc Rest[Rest[v]];
    v  = Join[{fdlo}, vi, {0.}];
    v  = MapThread[Max, {v, fdpay}];                 (* early exercise      *)
    k  = k + 1];
  v[[fdmid]]];

(* ---- 3. Monte-Carlo value at risk --------------------------------------- *)

vrn = 250000;      (* scenarios                                             *)
vrk = 64;          (* assets                                                *)
vrs = RandomReal[{-0.02, 0.02}, {vrn, vrk}];         (* timed: random       *)
vrw = Table[1./vrk, {vrk}];                          (* equal weights       *)
vrcs = Table[N[Mod[7 i + 3 j, 41] - 20]/1000., {i, 1, 20000}, {j, 1, vrk}];

(* Portfolio P&L is one dgemv; the risk number is a tail mean of the sorted
   result -- an order statistic, so this row is a Dot plus a Sort. *)
mcvar[s_, w_, q_] := Module[{pnl},
  pnl = Sort[s . w];
  Total[Take[pnl, q]]/q];

(* ---- run ---------------------------------------------------------------- *)

Print["Experiment 15 -- option pricing"];
Print[""];

bench["binomial American put, 4000 steps",        amtree[]];
check["binomial American put", amtree[]];

bench["explicit FD American put, 1000 x 25000",   amfd[]];
check["explicit FD American put", amfd[]];

bench["Monte-Carlo VaR, 250000 x 64",             mcvar[vrs, vrw, 10000]];
check["Monte-Carlo VaR", mcvar[vrcs, vrw, 200]];

Print[""];
Print["the two pricers agree to about 0.6%, which is what says they are"];
Print["pricing the same option rather than each reproducing itself"];

(* ---- the three spellings of the positive part --------------------------- *)
(*
   Until this experiment, Ramp did not exist in Mathilda and
   Clip[x, {0., Infinity}] returned UNEVALUATED because the infinite bound
   failed to parse -- so the only working spelling of max(x, 0) over an array
   was the two-pass mixed-type product.  All three work now; the row exists to
   show what the idiomatic one is worth.
*)
Print[""];
Print["-- the positive part, three spellings, 10^6 float64 --"];

pv = RandomReal[{-1, 1}, 10^6];
bench["Ramp[v]",                         Ramp[pv]];
bench["Clip[v, {0., Infinity}]",         Clip[pv, {0., Infinity}]];
bench["v UnitStep[v]",                   pv UnitStep[pv]];
