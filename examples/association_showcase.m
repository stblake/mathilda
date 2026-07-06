(* association_showcase.m — a realistic split-apply-combine pipeline over a
 * synthetic transaction log, using only the Association toolchain.
 *
 * Run:  ./Mathilda < examples/association_showcase.m
 * (companion write-up: examples/association-showcase.md)
 *
 * Every heavy step (Counts, GroupBy + reducer, ReverseSort, TakeLargest) is
 * hash-backed and amortised O(n), so the whole pipeline scales linearly. *)

SeedRandom[42];
categories = {"food", "rent", "travel", "fun", "misc"};
n = 100000;
txns = Table[{categories[[RandomInteger[{1, 5}]]], RandomInteger[{1, 100}]}, {n}];
Print["transactions: ", Length[txns]];

(* 1. How many transactions per category?  (hash tally, O(n)) *)
Print["counts:  ", Counts[txns[[All, 1]]]];

(* 2. Total spend per category  (split-apply-combine: group by category,
 *    reduce each group by summing its amounts). *)
spend = GroupBy[txns, First -> Last, Total];
Print["spend:   ", spend];

(* 3. Rank categories by total spend, descending. *)
Print["ranked:  ", ReverseSort[spend]];

(* 4. The single biggest-spending category. *)
Print["top-1:   ", TakeLargest[spend, 1]];

(* 5. Average transaction size per category (another reducer). *)
Print["avg:     ", GroupBy[txns, First -> Last, Mean] // N];

(* 6. Timings — confirm the pipeline is dominated by O(n) hashing. *)
Print["--- timings (ms) ---"];
Print["Counts:            ", 1000. Timing[Counts[txns[[All, 1]]]][[1]]];
Print["GroupBy+Total:     ", 1000. Timing[GroupBy[txns, First -> Last, Total]][[1]]];
Print["ReverseSort spend: ", 1000. Timing[ReverseSort[spend]][[1]]];
