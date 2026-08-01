(* ==========================================================================
   Experiment 13 -- Molecular dynamics: Lennard-Jones with a cut-off
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md for the
   measurements and the analysis.

       Mathilda      ./Mathilda -file molecular_dynamics.m
       Mathematica   wolframscript -file molecular_dynamics.m
       Python        python3 molecular_dynamics.py

   WHAT IT MEASURES.  The N-body row of the third sweep is the closest
   existing benchmark and it is a DIFFERENT kernel: gravity has no cut-off, so
   an all-pairs gravity step is pure arithmetic over a dense matrix.
   Lennard-Jones multiplies through a DATA-DEPENDENT mask,

       mk = UnitStep[rc2 - r2]

   so the question is whether a comparison result can stay on the machine
   buffer, or whether the mask forces the whole 2048 x 2048 interaction matrix
   back through the evaluator.

   DETERMINISM.  The configuration is a deterministically perturbed
   simple-cubic lattice, so the WHOLE benchmark reproduces across systems
   rather than only its check -- and the check is the potential energy, which
   every one of the 2 million pairs contributes to.

   THE SELF-INTERACTION.  A large constant is added to the diagonal of r2 so
   that i == j contributes nothing: 1/r2 is then 1e-12 rather than infinite,
   and the cut-off mask is 0 there in any case.  Writing it this way keeps the
   kernel branch-free, which is the whole point of the array formulation.
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

(* ---- configuration ------------------------------------------------------ *)

nmd   = 2048;        (* atoms                                               *)
mdrc2 = 6.25;        (* squared cut-off radius, 2.5 sigma                    *)
mddt  = 0.0005;      (* Verlet time step                                     *)

mdx = Table[N[Mod[i, 13]]                + 0.13 Sin[1.7 i], {i, 0, nmd - 1}];
mdy = Table[N[Mod[Quotient[i, 13], 13]]  + 0.13 Cos[2.3 i], {i, 0, nmd - 1}];
mdz = Table[N[Quotient[i, 169]]          + 0.13 Sin[0.9 i], {i, 0, nmd - 1}];

mdbig  = 1.0*^12 IdentityMatrix[nmd];    (* kills the self-interaction       *)
mdzero = Table[0., {nmd}];

(* ---- kernels ------------------------------------------------------------ *)

(* Lennard-Jones force on every atom.  Returns {fx, fy, fz}.

   NOTE that this does NOT share its intermediates with lje below, even though
   the two need the same six matrices.  Factoring them into a helper that
   returns {dx, dy, dz, i2, i6, mk} is the obvious way to write it and it is a
   TRAP in Mathilda: five of the six are float64 and the mask is int64, a
   mixed-dtype list of packed rows cannot be absorbed into one buffer, and the
   no-nesting invariant then materialises all six.  Measured below, and the
   subject of the README's main finding. *)
ljf[xs_, ys_, zs_] := Module[{dx, dy, dz, r2, mk, i2, i6, ff},
  dx = Outer[Subtract, xs, xs];
  dy = Outer[Subtract, ys, ys];
  dz = Outer[Subtract, zs, zs];
  r2 = dx dx + dy dy + dz dz + mdbig;
  mk = UnitStep[mdrc2 - r2];               (* 1 inside the cut-off, else 0   *)
  i2 = 1./r2;
  i6 = i2 i2 i2;                           (* (sigma/r)^6                    *)
  ff = (48. i6 i6 - 24. i6) i2 mk;         (* -dU/dr / r, masked             *)
  {Total[ff dx, {2}], Total[ff dy, {2}], Total[ff dz, {2}]}];

(* Total potential energy, 4 (r^-12 - r^-6), halved because every pair is
   counted twice.  This is the cross-system check: it is sensitive to every
   pair distance, so two systems agreeing on it are running the same physics. *)
lje[xs_, ys_, zs_] := Module[{dx, dy, dz, r2, mk, i2, i6},
  dx = Outer[Subtract, xs, xs];
  dy = Outer[Subtract, ys, ys];
  dz = Outer[Subtract, zs, zs];
  r2 = dx dx + dy dy + dz dz + mdbig;
  mk = UnitStep[mdrc2 - r2];
  i2 = 1./r2;
  i6 = i2 i2 i2;
  Total[4.(i6 i6 - i6) mk, 2]/2.];

(* One velocity-Verlet step.  State is {positions, velocities, forces}, each a
   3-element list of nmd-vectors -- uniform float64 throughout, so the whole
   state packs into one buffer per component. *)
mdstep[st_] := Module[{p, v, f, vh, pn, fn},
  p = st[[1]]; v = st[[2]]; f = st[[3]];
  vh = v + (0.5 mddt) f;                   (* half-kick                      *)
  pn = p + mddt vh;                        (* drift                          *)
  fn = ljf[pn[[1]], pn[[2]], pn[[3]]];
  {pn, vh + (0.5 mddt) fn, fn}];           (* second half-kick               *)

mdrun[m_] := Module[{st},
  st = Nest[mdstep,
            {{mdx, mdy, mdz}, {mdzero, mdzero, mdzero}, ljf[mdx, mdy, mdz]},
            m];
  lje[st[[1, 1]], st[[1, 2]], st[[1, 3]]]];

(* Cell-list binning: the irregular half of a real MD code.  100000 atoms are
   assigned to one of 8^3 cells, sorted by cell, and counted, which is what
   builds the neighbour-list offsets. *)
nmc = 100000;
mcx = Table[8. Mod[1.7 i, 1.], {i, 1, nmc}];
mcy = Table[8. Mod[2.3 i, 1.], {i, 1, nmc}];
mcz = Table[8. Mod[3.1 i, 1.], {i, 1, nmc}];

mdcells[] := Module[{ci, tl},
  ci = Floor[mcx] + 8 Floor[mcy] + 64 Floor[mcz] + 1;
  tl = Tally[Sort[ci]];
  Total[Accumulate[tl[[All, 2]]]]];

(* ---- run ---------------------------------------------------------------- *)

Print["Experiment 13 -- molecular dynamics"];
Print[""];

bench["Lennard-Jones force, 2048 atoms, cut-off", ljf[mdx, mdy, mdz]];
check["Lennard-Jones force", lje[mdx, mdy, mdz]];

bench["velocity-Verlet, 2048 atoms, 10 steps",    mdrun[10]];
check["velocity-Verlet", mdrun[10]];

bench["cell-list binning, 100000 atoms",          mdcells[]];
check["cell-list binning", mdcells[]];

(* ---- the finding: a mixed-dtype tuple return ---------------------------- *)
(*
   Returning several arrays from one function is free when they share an
   element type and ruinous when they do not.  In Mathilda a list of packed
   rows is absorbed into one buffer only if the rows agree in dtype; a float64
   / int64 mixture declines, and the no-nesting invariant then materialises
   EVERY element.  The caller pays on its next operation, not at the return,
   which is what makes it hard to find by profiling.

   In Mathematica both forms are cheap, so this section is a Mathilda finding
   that the Mathematica column exists to prove is not intrinsic.
*)
Print[""];
Print["-- the mixed-dtype tuple return (600 x 600) --"];

tta = RandomReal[{0, 1}, {600, 600}];
ttb = RandomReal[{0, 1}, {600, 600}];
ttm = UnitStep[tta - 0.5];               (* int64 mask                       *)
ttf = 1. UnitStep[tta - 0.5];            (* the same mask, widened to Real   *)

bench["{a, b}          -- uniform float64",       {tta, ttb}];
bench["{a, b, mask}    -- one int64 row",         {tta, ttb, ttm}];
bench["{a, b, 1. mask} -- widened by hand",       {tta, ttb, ttf}];

ttu = {tta, ttb};
ttx = {tta, ttb, ttm};
bench["caller's next op after the uniform tuple", ttu[[1]] + ttu[[2]]];
bench["caller's next op after the mixed tuple",   ttx[[1]] + ttx[[2]]];
