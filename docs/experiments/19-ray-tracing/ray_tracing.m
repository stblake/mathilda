(* ==========================================================================
   Experiment 19 -- Ray tracing: branch-free arrays, and a 65-element table
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.

       Mathilda      ./Mathilda -file ray_tracing.m
       Mathematica   wolframscript -file ray_tracing.m
       Python        python3 ray_tracing.py

   WHAT IT MEASURES.  A ray tracer written for an array language is
   BRANCH-FREE by construction: there is no If anywhere, because every
   conditional becomes a mask and every selection becomes arithmetic on that
   mask.  That style is how all high-throughput array code expresses control
   flow, and the suite had never measured it end to end.

   It also ends with the operation experiment 12 is about: once the winning
   sphere per ray is known, its centre must be GATHERED by index -- from a
   65-entry table, with 262144 indices.

   DETERMINISM.  The scene, the camera and the light are all closed forms, so
   the mean pixel intensity is an exact cross-system check.
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

(* ---- camera and scene --------------------------------------------------- *)

rtw = 512; rth = 512;
rtnp = rtw rth;                (* 262144 rays, one per pixel                *)
rtns = 64;                     (* spheres                                   *)
rtr  = 0.35; rtr2 = rtr^2;

(* Every ray is flattened into a plain vector, so the whole image is one
   array operation per sphere rather than a loop over pixels. *)
rtjs = Flatten[Table[(j - (rtw + 1)/2.)/rtw, {i, 1, rth}, {j, 1, rtw}]];
rtis = Flatten[Table[(i - (rth + 1)/2.)/rth, {i, 1, rth}, {j, 1, rtw}]];
rtnm = 1./Sqrt[rtjs^2 + rtis^2 + 1.];              (* normalise             *)
rtdx = rtjs rtnm; rtdy = rtis rtnm; rtdz = rtnm;

rtcx = Table[-3.5 + Mod[s - 1, 8],      {s, 1, rtns}];
rtcy = Table[-3.5 + Quotient[s - 1, 8], {s, 1, rtns}];
rtcz = Table[6. + 0.4 Sin[1.7 s],       {s, 1, rtns}];

(* A dummy entry at position 1 so a MISS (sphere id 0) can index the table
   without a branch: cid + 1 is then always in range, and the miss is masked
   out at the end instead. *)
rtcxp = Prepend[rtcx, 0.];
rtcyp = Prepend[rtcy, 0.];
rtczp = Prepend[rtcz, 0.];

rtlx = 0.577; rtly = -0.577; rtlz = -0.577;        (* light direction       *)

(* ---- the tracer --------------------------------------------------------- *)

rtrender[] := Module[
  {tmin, cid, s, cx, cy, cz, tca, d2, disc, tt, ok, cand, bt,
   hx, hy, hz, gx, gy, gz, nx, ny, nz, lam, vis},

  tmin = Table[1.0*^9, {rtnp}];          (* nearest hit so far, per ray     *)
  cid  = Table[0, {rtnp}];               (* which sphere; 0 means a miss    *)
  s    = 1;

  While[s <= rtns,
    cx = rtcx[[s]]; cy = rtcy[[s]]; cz = rtcz[[s]];

    (* Ray-sphere intersection, from the eye at the origin.  tca is the
       projection of the centre onto the ray; disc < 0 means a miss. *)
    tca  = rtdx cx + rtdy cy + rtdz cz;
    d2   = (cx cx + cy cy + cz cz) - tca tca;
    disc = rtr2 - d2;

    (* Clip before Sqrt so a miss produces 0 rather than a NaN: the mask
       below discards it, but a NaN would propagate through the arithmetic. *)
    tt = tca - Sqrt[Clip[disc, {0., 1.0*^9}]];

    ok = UnitStep[disc] UnitStep[tt - 0.001];      (* hit, and in front     *)

    (* A miss is pushed to +infinity so the running minimum ignores it. *)
    cand = ok tt + (1 - ok) 1.0*^9;

    (* STRICTLY closer: the epsilon matters, because at s = 1 a missing ray
       has cand == tmin == 1e9 exactly, and UnitStep[0] is 1 -- so without it
       every missed ray would claim to have hit sphere 1. *)
    bt   = UnitStep[tmin - cand - 1.0*^-6];
    cid  = bt s + (1 - bt) cid;                    (* select, no branch     *)
    tmin = MapThread[Min, {tmin, cand}];
    s    = s + 1];

  hx = rtdx tmin; hy = rtdy tmin; hz = rtdz tmin;  (* hit points            *)

  (* THE GATHER: 262144 indices into a 65-entry table. *)
  gx = rtcxp[[cid + 1]]; gy = rtcyp[[cid + 1]]; gz = rtczp[[cid + 1]];

  nx = (hx - gx)/rtr; ny = (hy - gy)/rtr; nz = (hz - gz)/rtr;   (* normals  *)
  lam = Clip[nx rtlx + ny rtly + nz rtlz, {0., 1.}];            (* Lambert  *)
  vis = UnitStep[1.0*^8 - tmin];                                (* hit mask *)
  Total[lam vis]/rtnp];

(* ---- run ---------------------------------------------------------------- *)

Print["Experiment 19 -- ray tracing"];
Print[""];

bench["ray trace, 512^2 rays x 64 spheres, diffuse", rtrender[]];
check["ray trace (mean pixel intensity)", rtrender[]];

(* ---- the gather, alone -------------------------------------------------- *)
(*
   The source table has 65 entries -- small BY DESIGN, not by accident, so no
   packing threshold will ever pack it -- and the index array has 262144.
   Until Part learned to lift a below-threshold source, this was one boxed
   expression per output element.
*)
Print[""];
Print["-- the sphere-parameter lookup, alone --"];

rtcid = RandomInteger[{0, 64}, rtnp];
bench["table[[cid + 1]], 262144 indices from 65", rtcxp[[rtcid + 1]]];

(* One sphere's worth of the inner loop, to show where the rest of the time
   goes: fifteen whole-array passes over 262144 float64 is ~32 MB of traffic
   per sphere, and there are 64 spheres. *)
Print[""];
Print["-- one sphere of the inner loop --"];
rttca = rtdx 1. + rtdy 2. + rtdz 6.;
bench["tca = dx cx + dy cy + dz cz",  rtdx 1. + rtdy 2. + rtdz 6.];
bench["Sqrt[Clip[disc, ...]]",        Sqrt[Clip[rtr2 - rttca, {0., 1.0*^9}]]];
bench["UnitStep[disc] UnitStep[...]", UnitStep[rttca] UnitStep[rttca - 0.001]];
bench["MapThread[Min, {a, b}]",       MapThread[Min, {rttca, rtdz}]];
