(* ==========================================================================
   Experiment 12 -- Graph analytics and sparse linear algebra
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md for the
   measurements and the analysis; this file is the code those measurements
   were taken from.

       Mathilda      ./Mathilda -file graph_and_sparse.m
       Mathematica   wolframscript -file graph_and_sparse.m
       Python        python3 graph_and_sparse.py

   WHAT IT MEASURES.  Every sparse and graph kernel is built out of one
   operation: gather -- y = x[[idx]], reading a value array through an index
   array.  A sparse matrix-vector product is a gather, a multiply and a
   segmented sum; a PageRank iteration is twenty of those; a breadth-first
   search is a gather of adjacency rows and a set difference.

   REPRESENTATION.  Mathilda has no SparseArray, so the graph is stored the
   way a system without one actually stores it, and the way every
   uniform-degree graph arrives anyway: a dense gn x gdeg matrix of neighbour
   ids.  That is CSR with the row pointers implied, and it isolates the
   gather rather than hiding it inside a sparse-matrix type.

   DETERMINISM.  The timed graph is random -- 100000 nodes, degree 16, so
   1.6e6 edges, the size of a mid-range citation or road network.  The three
   systems cannot be made to draw the same random graph, so every value check
   runs a separate DETERMINISTIC 4096-node graph instead.  A timing is
   meaningless until the three systems agree on the check.
   ========================================================================== *)

(* ---- shared reporting helpers (identical in every experiment file) ------- *)

SetAttributes[bench, HoldRest];

(* bench[label, expr] -- one untimed warm-up, then the MINIMUM of three timed
   runs.  The minimum, not the mean: we are measuring the cost of the work,
   and every source of noise on a loaded machine can only add. *)
bench[label_String, expr_] := Module[{ts},
  expr;
  ts = Table[First[AbsoluteTiming[expr]], {3}];
  Print[StringPadRight[label, 52], ToString[Round[1000. Min[ts], 0.001]], " ms"]
];

check[label_String, value_] :=
  Print[StringPadRight[label, 52], "check = ", value];

(* ---- the timed graph ---------------------------------------------------- *)

gn   = 100000;                        (* nodes                               *)
gdeg = 16;                            (* in-degree, uniform                  *)

(* gadj[[i]] is the list of nodes that link TO node i.  Uniform in-degree is
   what makes the row pointers implicit: row i occupies the flat positions
   (i-1) gdeg + 1 .. i gdeg. *)
gadj  = RandomInteger[{1, gn}, {gn, gdeg}];
gflat = Flatten[gadj];                (* the CSR column-index array          *)
gv    = RandomReal[{0, 1}, gn];       (* a value per node                    *)
gw    = RandomReal[{0, 1}, gn gdeg];  (* a weight per edge                   *)

(* ---- the deterministic check graph -------------------------------------- *)

gcn = 4096;
gcd = 8;
gca = Table[Mod[i j + j, gcn] + 1, {i, 1, gcn}, {j, 1, gcd}];
gcf = Flatten[gca];
gcw = Table[N[Mod[i, 7] + 1]/7., {i, 1, gcn gcd}];

(* ---- kernels ------------------------------------------------------------ *)

(* PageRank by power iteration.
     p[i] <- (1-d)/n + d * sum over the in-neighbours j of i of p[j]/deg

   Written as three whole-array operations, which is the point: the gather
   p[[f]] produces one value per EDGE, Partition reshapes those into one row
   per node, and Total[..., {2}] sums each row.  There is no scalar loop and
   no sparse-matrix type -- just the gather and a segmented reduction. *)
pgrk[f_, nn_, dd_, m_] := Module[{p, k},
  p = Table[1./nn, {nn}];
  k = 0;
  While[k < m,
    p = 0.15/nn + 0.85 Total[Partition[p[[f]], dd], {2}]/dd;
    k = k + 1];
  p];

(* Breadth-first search: level-synchronous frontier expansion.

   a[[fr]] gathers whole adjacency ROWS -- a rank-2 gather -- and the set
   operations then do the work of a visited-bitmap.  This is the Graph500
   kernel written the way an array language expresses it. *)
bfsr[a_, st_, lev_] := Module[{vis, fr, nx, k},
  vis = {st}; fr = {st}; k = 0;
  While[k < lev && Length[fr] > 0,
    nx  = Complement[Union[Flatten[a[[fr]]]], vis];
    vis = Union[vis, nx];
    fr  = nx;
    k   = k + 1];
  Length[vis]];

(* ---- run ---------------------------------------------------------------- *)

Print["Experiment 12 -- graph analytics and sparse linear algebra"];
Print[""];

(* 1. The primitive, alone: 1.6e6 indices into a 100000-element vector. *)
bench["gather, 1.6e6 indices into 100000", gv[[gflat]]];
check["gather", Total[gcw[[Table[Mod[7 i, 4096] + 1, {i, 1, 32768}]]]]];

(* 2. One sparse matrix-vector product: gather, multiply, segmented sum. *)
bench["SpMV, 100000 x 16 CSR", Total[Partition[gw gv[[gflat]], gdeg], {2}]];
check["SpMV", Total[Total[Partition[gcw gcw[[gcf]], gcd], {2}]]];

(* 3. Twenty of the above, with the damping term. *)
bench["PageRank, 1.6e6 edges, 20 iterations", pgrk[gflat, gn, gdeg, 20]];
check["PageRank", Total[pgrk[gcf, gcn, gcd, 20]]];

(* 4. Five levels of frontier expansion from one source. *)
bench["breadth-first search, 5 levels", bfsr[gadj, 1, 5]];
check["breadth-first search", bfsr[gca, 1, 5]];
