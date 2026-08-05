(* Experiment 23 -- The per-operation constant.
   ROADMAP ITEM 7 -- "~8 us per array operation, which dominates any loop with a
   short body", named as "the whole gap on three experiments" (15, 16, 18).

   METHOD.  Hold the total element count fixed and vary the OPERATION count.  A
   system with zero per-op overhead gives the same time for 1 op on 10^6 elements
   as for 1000 ops on 10^3.  The gap between those two rows IS the constant, and
   dividing it by the op count gives it in microseconds directly.

   This is the one experiment whose finding is a NUMBER rather than a ratio, and
   it is why it deserves its own file rather than a note in another. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Total", "Do", "Table", "Plus", "Times", "Length"}];

big = rand01[{1000000}];
small = rand01[{1000}];
tiny = rand01[{10}];

(* 1 operation, 10^6 elements: overhead is negligible, this is bandwidth. *)
bench["1 op on 10^6 elements", Total[big];];
check["1 op on 10^6 elements", Round[10^6 N[Total[Range[10]]]]];

(* 1000 operations, 10^3 elements each: same element count, 1000x the overhead. *)
bench["1000 ops on 10^3 elements", Do[Total[small], {1000}];];
check["1000 ops on 10^3 elements", Round[10^6 N[Total[Range[10]]]]];

(* 100000 operations on 10 elements: overhead completely dominates. *)
bench["100000 ops on 10 elements", Do[Total[tiny], {100000}];];
check["100000 ops on 10 elements", Round[10^6 N[Total[Range[10]]]]];

(* The same sweep for an elementwise binary op rather than a reduction. *)
bench["1 elementwise op on 10^6", (big + big);];
check["1 elementwise op on 10^6", Round[10^6 N[Total[Range[10] + Range[10]]]]];

bench["10000 elementwise ops on 10^2",
  Module[{s = rand01[{100}]}, Do[s + s, {10000}]];];
check["10000 elementwise ops on 10^2",
  Round[10^6 N[Total[Range[10] + Range[10]]]]];

(* A scalar loop: the pure interpreter dispatch cost, no arrays at all. *)
bench["100000 scalar additions", Module[{s = 0.}, Do[s = s + 1., {100000}]];];
check["100000 scalar additions", 1];
