(* ==========================================================================
   Experiment 14 -- Bioinformatics: the integer half of the buffer
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.

       Mathilda      ./Mathilda -file sequence_alignment.m
       Mathematica   wolframscript -file sequence_alignment.m
       Python        python3 sequence_alignment.py

   WHAT IT MEASURES.  Everything here is INTEGER end to end.  Bases are
   symbols, scores are small integers, gap penalties are integers, k-mers are
   integers, and an alignment matrix is integer dynamic programming.  If the
   int64 arms of the elementwise, scan and structural paths are missing, an
   integer workload finds it and a floating-point one never will.

   THE ALIGNMENT, AND WHY IT IS WRITTEN THIS WAY.  The Needleman-Wunsch
   within-row gap recurrence is

       h[j] = max(cand[j], h[j-1] - g)

   which is a max-plus prefix scan and looks sequential.  Substituting
   H[j] = h[j] + g j turns it into a plain running maximum:

       H[j] = max(cand[j] + g j, H[j-1])          h[j] = H[j] - g j

   so all three systems run THAT -- FoldList[Max, ...] in the two CAS,
   np.maximum.accumulate in NumPy.  The columns then compare execution rather
   than cleverness.  This is the standard vectorisation and is what a
   performance-minded implementer writes.

   DETERMINISM.  Both sequences and the k-mer input are generated from closed
   forms, so every value is an integer that all three systems must agree on
   exactly -- no tolerance is involved.
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

(* ---- 1. Needleman-Wunsch global alignment ------------------------------- *)

bm   = 2000;      (* sequence length; the DP matrix is bm x bm              *)
bgap = 1;         (* linear gap penalty                                     *)

bs1 = Table[Mod[i^2 + 3 i, 4], {i, 1, bm}];              (* bases 0..3       *)
bs2 = Table[Mod[7 i + Quotient[i, 3], 4], {i, 1, bm}];

brng1 = Range[bm];         (* 1..bm, the g j shift                          *)
brng0 = Range[0, bm];      (* 0..bm, the g j unshift                        *)

(* One row of the DP matrix, from the row above it.

   `sco` is the substitution score, +2 on a match and -1 on a mismatch,
   computed WITHOUT a comparison predicate: for integer bases, Abs[b - ch] is
   0 exactly when they are equal, so UnitStep[Abs[b - ch] - 1] is the
   "not equal" indicator and 2 - 3 (not equal) is the score.  Written with
   Map[If[...]] instead it would be one interpreted call per column. *)
nwrow[prev_, ch_] := Module[{sco, dg, up, cand, hh},
  sco  = 2 - 3 UnitStep[Abs[bs2 - ch] - 1];
  dg   = Most[prev] + sco;                    (* diagonal: match/mismatch    *)
  up   = Rest[prev] - bgap;                   (* from above: a gap           *)
  cand = MapThread[Max, {dg, up}] + bgap brng1;
  hh   = FoldList[Max, prev[[1]] - bgap, cand];   (* the running maximum     *)
  hh - bgap brng0];                               (* undo the shift          *)

nwalign[] := Last[Fold[nwrow, -bgap brng0, bs1]];

(* ---- 2. k-mer encoding and distinct count ------------------------------- *)

bkn = 500000;     (* bases                                                   *)
bk  = 12;         (* k-mer length; 4^12 = 16.7e6 possible keys               *)

bcode = Table[Mod[i^2 + 3 i + Quotient[i, 7], 4], {i, 1, bkn}];
bpow  = Table[4^(bk - j), {j, 1, bk}];       (* base-4 place values          *)

(* Partition[..., bk, 1] is the SLIDING WINDOW -- offset 1, so the rows
   overlap.  The Dot that follows is the base-4 encoding of each window into a
   single integer key, and is a matrix-vector product over exact integers. *)
kmers[] := Length[Union[Partition[bcode, bk, 1] . bpow]];

(* ---- 3. rolling GC content ---------------------------------------------- *)

bgn  = 10^7;
bgw  = 1000;                                  (* window                      *)
bgseq = RandomInteger[{0, 1}, bgn];           (* timed input: random         *)
bgcs  = Table[Mod[i^3 + i + Quotient[i, 3], 2], {i, 1, 100000}];  (* check    *)

(* A rolling sum over a window is two shifted reads of the prefix sum -- the
   standard trick, and the reason this row is a scan benchmark rather than a
   convolution one. *)
gcwin[s_, w_] := Module[{cs},
  cs = Accumulate[s];
  Total[Drop[cs, w] - Drop[cs, -w]]];

(* ---- run ---------------------------------------------------------------- *)

Print["Experiment 14 -- sequence alignment and k-mer counting"];
Print[""];

bench["Needleman-Wunsch alignment, 2000 x 2000", nwalign[]];
check["Needleman-Wunsch alignment", nwalign[]];

bench["k-mer encode + distinct count, 5e5, k=12", kmers[]];
check["k-mer encode + distinct count", kmers[]];

bench["rolling GC content, 10^7 bases, w=1000",   gcwin[bgseq, bgw]];
check["rolling GC content", gcwin[bgcs, 100]];

(* ---- the primitives underneath ------------------------------------------ *)
(*
   The alignment row is built out of five operations on 2000-element int64
   vectors.  Timing them separately is what found this experiment's result:
   two of them had no integer arm at all, and were handing unpacked operands
   to the three that did.  Run this section against an older Mathilda to see
   Abs and Most at 460 and 138 microseconds.
*)
Print[""];
Print["-- the inner row's primitives, per 2000-element call --"];

nwprev = Table[-i, {i, 0, bm}];

bench["Abs[seq - ch]",                    Do[Abs[bs2 - 2], {200}]];
bench["UnitStep[...]",                    Do[UnitStep[bs2 - 2], {200}]];
bench["Most[prev]",                       Do[Most[nwprev], {200}]];
bench["MapThread[Max, {Most, Rest}]",     Do[MapThread[Max, {Most[nwprev], Rest[nwprev]}], {200}]];
bench["FoldList[Max, 0, v]",              Do[FoldList[Max, 0, brng1], {200}]];
bench["the whole row",                    Do[nwrow[nwprev, 2], {200}]];

Print[""];
Print["-- the k-mer pipeline, split --"];
bench["Partition[code, 12, 1] (the window)", Partition[bcode, bk, 1]];
kmpm = Partition[bcode, bk, 1];
bench["... . pow4 (the encode)",             kmpm . bpow];
kmkeys = kmpm . bpow;
bench["Union (the distinct count)",          Union[kmkeys]];
