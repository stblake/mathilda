(* Experiment 28 -- The start/step/n selector.
   ROADMAP ITEM 6 -- "a start/step/n selector: no position array for spans,
   strided writes or gathers", valued at 2x on spans, most of the graph rows, and
   the sieve (experiments 5, 11, 12).

   THE DEFECT.  A strided span v[[1;;n;;2]] is describable in three integers, but
   if the implementation materialises a position array first it pays n boxed
   integers before it reads a single element.  The sieve is the extreme case: it
   is nothing BUT strided writes, and experiment 9 measured it 25x behind NumPy
   while being 1.18x ahead of Mathematica -- the two baselines disagreeing is
   what identified this as a machine-level gap rather than a competitive one. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Span", "Part", "Range", "Table", "ReplacePart", "Take", "Drop"}];

n = 4000000;
v = rand01[{n}];

bench["span v[[1;;n;;2]]", v[[1 ;; n ;; 2]];];
check["span v[[1;;n;;2]]", Total[Range[10][[1 ;; 10 ;; 2]]]];

bench["span v[[1;;n;;7]]", v[[1 ;; n ;; 7]];];
check["span v[[1;;n;;7]]", Total[Range[20][[1 ;; 20 ;; 7]]]];

bench["span v[[2;;n;;2]] (offset)", v[[2 ;; n ;; 2]];];
check["span v[[2;;n;;2]] (offset)", Total[Range[10][[2 ;; 10 ;; 2]]]];

bench["contiguous span v[[1;;n/2]]", v[[1 ;; 2000000]];];
check["contiguous span v[[1;;n/2]]", Total[Range[10][[1 ;; 5]]]];

(* Range itself: generating n integers with a step, no source array at all. *)
bench["Range[1, 4x10^6, 3]", Range[1, 4000000, 3];];
check["Range[1, 4x10^6, 3]", Total[Range[1, 20, 3]]];

(* The sieve of Eratosthenes: nothing but strided writes. *)
sieve[lim_] := Module[{s = ConstantArray[1, lim], i = 2},
  While[i i <= lim,
    If[s[[i]] === 1,
      s = ReplacePart[s, Table[{j} -> 0, {j, i i, lim, i}]]];
    i++];
  Total[s] - 1];
bench["sieve to 200000 (strided writes)", sieve[200000], 1];
check["sieve to 200000 (strided writes)", sieve[100]];
