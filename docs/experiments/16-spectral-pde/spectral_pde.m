(* ==========================================================================
   Experiment 16 -- Spectral PDE: an FFT inside a time loop
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.

       Mathilda      ./Mathilda -file spectral_pde.m
       Mathematica   wolframscript -file spectral_pde.m
       Python        python3 spectral_pde.py

   WHAT IT MEASURES.  The suite's existing spectral rows time ONE large
   transform, so what they measure is FFTW.  A pseudo-spectral PDE solver has
   the opposite profile: thousands of MEDIUM transforms with complex
   elementwise algebra between them.  What that measures is per-call overhead,
   complex-array arithmetic, and whether a Fourier result stays on the machine
   buffer for the next operation.

   THE TRANSFORM CONVENTION IS LOAD-BEARING.  Wolfram's Fourier is
   (1/Sqrt[n]) Sum x e^{+2 Pi i ...} -- a POSITIVE exponent -- which is NumPy's
   ifft(norm="ortho"); InverseFourier is fft(norm="ortho").  Getting it
   backwards gives a plausible-looking wrong answer: the solution still
   evolves, still looks like turbulence, and is a different equation.  The
   cross-system value check is what stops that, and it is worth more here than
   anywhere else in the suite.

   Because u = InverseFourier[uh] carries e^{-2 Pi i m x / L}, differentiation
   in coefficient space is multiplication by -I kappa with kappa = 2 Pi m / L.

   THE CHECKS ARE SHORT RUNS ON PURPOSE.  Both equations are chaotic, so the
   state at the final time is a property of the arithmetic rather than of the
   equation; asking three independently written solvers to agree on it to six
   figures would be meaningless.  50 KS steps and 20 NS steps are still well
   inside the deterministic regime.
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

(* ---- 1. Kuramoto-Sivashinsky, u_t = -u u_x - u_xx - u_xxxx -------------- *)

ksn  = 2048;
ksll = 64. Pi;         (* domain length; large enough to be chaotic          *)
ksdt = 0.01;

ksxx = Table[ksll (j - 1)/ksn, {j, 1, ksn}];

(* Wavenumbers in FFT order: 0, 1, ..., n/2-1, -n/2, ..., -1, scaled to the
   physical domain. *)
kskk = Table[(2. Pi/ksll) If[j <= ksn/2, j - 1, j - 1 - ksn], {j, 1, ksn}];

(* Linear operator: -u_xx - u_xxxx becomes kappa^2 - kappa^4.  It is stiff
   (kappa^4 reaches 2.7e8), which is why the scheme below treats it
   implicitly -- an explicit one would need an unusably small step. *)
kslin = kskk^2 - kskk^4;
ksden = 1./(1. - ksdt kslin);

ksu0 = Cos[ksxx/16.](1. + Sin[ksxx/16.]);
ksh0 = Fourier[ksu0];

(* IMEX Euler: the stiff linear part implicit, the nonlinear part explicit.
   -u u_x is written as -1/2 d/dx (u^2), which costs one transform instead of
   two and is exactly conservative. *)
ksstep[uh_] := Module[{u, nl},
  u  = Re[InverseFourier[uh]];
  nl = (0.5 I kskk) Fourier[u u];
  (uh + ksdt nl) ksden];

ksrun[m_] := Total[Abs[Nest[ksstep, ksh0, m]]^2];

(* ---- 2. 2D Navier-Stokes in vorticity form ------------------------------ *)

nsn  = 128;
nsdt = 0.002;
nsnu = 0.002;          (* viscosity                                          *)

nsk1  = Table[If[j <= nsn/2, j - 1, j - 1 - nsn], {j, 1, nsn}];
nskx  = Table[nsk1[[j]], {i, 1, nsn}, {j, 1, nsn}];
nsky  = Table[nsk1[[i]], {i, 1, nsn}, {j, 1, nsn}];
nsk2  = nskx^2 + nsky^2;

(* The (0,0) mode is the mean, which the Poisson solve cannot determine; add 1
   there so the reciprocal is finite, and the mode is never used. *)
nsk2i = 1./(nsk2 + Table[If[i == 1 && j == 1, 1., 0.], {i, 1, nsn}, {j, 1, nsn}]);

(* 2/3 dealiasing: a quadratic nonlinearity aliases the top third of the
   spectrum back onto the resolved modes, so those coefficients are zeroed. *)
nsdeal = Table[If[nsk1[[i]]^2 + nsk1[[j]]^2 <= (nsn/3)^2, 1., 0.],
               {i, 1, nsn}, {j, 1, nsn}];

nsw0 = Fourier[Table[Sin[2. Pi i/nsn] Cos[4. Pi j/nsn]
                     + 0.4 Sin[6. Pi (i + j)/nsn], {i, 1, nsn}, {j, 1, nsn}]];

(* Stream function from vorticity by a Poisson solve, velocity from the
   stream function, then the advective term -- six transforms per step. *)
nsstep[wh_] := Module[{ph, u, v, wx, wy, nl},
  ph = wh nsk2i;                                    (* -Laplacian^-1 w       *)
  u  = Re[InverseFourier[(-I nsky) ph]];            (*  d psi / dy           *)
  v  = Re[InverseFourier[( I nskx) ph]];            (* -d psi / dx           *)
  wx = Re[InverseFourier[(-I nskx) wh]];
  wy = Re[InverseFourier[(-I nsky) wh]];
  nl = Fourier[u wx + v wy] nsdeal;
  (wh - nsdt nl)/(1. + nsdt nsnu nsk2)];            (* implicit viscosity    *)

nsrun[m_] := Total[Abs[Nest[nsstep, nsw0, m]]^2, 2];

(* ---- 3. FFT Poisson solve ----------------------------------------------- *)

psn   = 512;
pssrc = Table[Sin[2. Pi i/psn] Cos[6. Pi j/psn], {i, 1, psn}, {j, 1, psn}];
psk1  = Table[If[j <= psn/2, j - 1, j - 1 - psn], {j, 1, psn}];
psk2  = Table[N[psk1[[i]]^2 + psk1[[j]]^2], {i, 1, psn}, {j, 1, psn}];
psk2i = 1./(psk2 + Table[If[i == 1 && j == 1, 1., 0.], {i, 1, psn}, {j, 1, psn}]);

(* The source is varied per iteration so that no cache can answer twice. *)
psolve[m_] := Module[{r, k},
  r = 0.; k = 0;
  While[k < m,
    r = Total[Abs[InverseFourier[Fourier[pssrc (1. + 0.001 k)] psk2i]]^2, 2];
    k = k + 1];
  r];

(* ---- run ---------------------------------------------------------------- *)

Print["Experiment 16 -- spectral PDE"];
Print[""];

bench["Kuramoto-Sivashinsky, 2048 modes, 2000 steps", ksrun[2000]];
check["Kuramoto-Sivashinsky (50 steps)", ksrun[50]];

bench["2D Navier-Stokes vorticity, 128^2, 200 steps", nsrun[200]];
check["2D Navier-Stokes (20 steps)", nsrun[20]];

bench["FFT Poisson solve, 512^2, 30 solves",          psolve[30]];
check["FFT Poisson solve (2 solves)", psolve[2]];

(* ---- where a time step's cost actually is -------------------------------- *)
(*
   The transforms are not the cost.  Two 2048-point transforms are about 30
   microseconds together and ten elementwise passes over 2048 complexes about
   4 microseconds of memory traffic; the rest of a ~166 microsecond step is
   evaluation overhead.  This section splits it so that claim is checkable.
*)
Print[""];
Print["-- one KS step, split --"];

kszh = ksh0;
bench["InverseFourier (one transform)", InverseFourier[kszh]];
kszu = Re[InverseFourier[kszh]];
bench["Fourier (one transform)",        Fourier[kszu kszu]];
bench["one complex elementwise pass",   kszh ksden];
bench["the whole step",                 ksstep[kszh]];
