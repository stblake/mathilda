#include "list_common.h"
#include "list.h"

void list_init(void) {
    symtab_add_builtin("Table", builtin_table);
    symtab_add_builtin("Range", builtin_range);
    symtab_add_builtin("Array", builtin_array);
    symtab_add_builtin("ConstantArray", builtin_constant_array);
    symtab_add_builtin("ArrayFlatten", builtin_array_flatten);
    symtab_add_builtin("Take", builtin_take);
    symtab_add_builtin("Drop", builtin_drop);
    symtab_add_builtin("Flatten", builtin_flatten);
    symtab_add_builtin("FlattenAt", builtin_flatten_at);
    symtab_add_builtin("Partition", builtin_partition);
    symtab_add_builtin("Pick", builtin_pick);
    symtab_get_def("Pick")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Pick",
        "Pick[expr, sel]\n\tPicks out the elements of expr for which the\n"
        "\tcorresponding element of sel is True.\n"
        "Pick[expr, sel, patt]\n\tPicks out the elements of expr for which the\n"
        "\tcorresponding element of sel matches patt.\n"
        "\tOperates at all levels; sel must mirror the structure of expr, and\n"
        "\tthe head of expr is preserved. Returns unevaluated if the structures\n"
        "\tdisagree.");
    symtab_add_builtin("RotateLeft", builtin_rotateleft);
    symtab_add_builtin("RotateRight", builtin_rotateright);
    symtab_add_builtin("Reverse", builtin_reverse);
    symtab_add_builtin("Rescale", builtin_rescale);
    symtab_add_builtin("PadRight", builtin_padright);
    symtab_add_builtin("PadLeft", builtin_padleft);
    symtab_add_builtin("Join", builtin_join);
    symtab_add_builtin("Catenate", builtin_catenate);
    symtab_get_def("Catenate")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Catenate",
        "Catenate[{e1, e2, ...}]\n\tConcatenates the ei (which must share a head)\n"
        "\tinto one, flattening a single level. A list of associations merges into\n"
        "\tone association (later keys win).");
    symtab_get_def("Join")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Join",
        "Join[list1, list2, ...]\n"
        "\tConcatenates lists or other expressions that share the same head.\n"
        "Join[list1, list2, ..., n]\n"
        "\tJoins the objects at level n in each of the lists.\n"
        "\tHandles ragged arrays by concatenating successive elements at level n.");
    symtab_add_builtin("Transpose", builtin_transpose);
    symtab_add_builtin("ConjugateTranspose", builtin_conjugate_transpose);
    symtab_add_builtin("Tally", builtin_tally);
    symtab_add_builtin("Union", builtin_union);
    symtab_add_builtin("Intersection", builtin_intersection);
    symtab_add_builtin("Complement", builtin_complement);
    symtab_add_builtin("DeleteDuplicates", builtin_deleteduplicates);
    symtab_add_builtin("DeleteDuplicatesBy", builtin_deleteduplicatesby);
    symtab_get_def("DeleteDuplicatesBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DeleteDuplicatesBy",
        "DeleteDuplicatesBy[expr, f]\n\tKeeps the first element for each distinct\n"
        "\tf[element], preserving order. Over an association, f is applied to the\n"
        "\tvalues and the surviving entries are kept (keys preserved).");
    symtab_add_builtin("Split", builtin_split);
    symtab_add_builtin("SplitBy", builtin_splitby);
    symtab_set_docstring("SplitBy",
        "SplitBy[list, f]\n"
        "\tsplits list into runs of consecutive elements that give the same\n"
        "\tvalue of f[element]. Only adjacent elements are grouped (unlike\n"
        "\tGatherBy, which collects equal keys from anywhere in the list).\n"
        "SplitBy[list, {f1, f2, ...}]\n"
        "\tsplits by f1, then splits each resulting run by f2, and so on,\n"
        "\tnesting one level deeper per function.");
    symtab_add_builtin("Gather", builtin_gather);
    symtab_get_def("Gather")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Gather",
        "Gather[list]\n"
        "\tGathers identical elements of list into sublists, giving\n"
        "\t{{group1}, {group2}, ...}. Sublists appear in order of the first\n"
        "\toccurrence of their element, and elements keep their input order\n"
        "\twithin a sublist. Equal elements are collected from anywhere in the\n"
        "\tlist, not only from adjacent runs (unlike Split).\n"
        "\tGather[list] is equivalent to GatherBy[list, Identity].");
    symtab_add_builtin("Total", builtin_total);
    symtab_add_builtin("Accumulate", builtin_accumulate);
    symtab_add_builtin("Differences", builtin_differences);
    symtab_add_builtin("Ratios", builtin_ratios);
    symtab_add_builtin("Commonest", builtin_commonest);
    symtab_add_builtin("Min", builtin_min);
    symtab_add_builtin("Max", builtin_max);
    symtab_add_builtin("MinMax", builtin_minmax);
    symtab_get_def("MinMax")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MinMax",
        "MinMax[list]\n\tGives {Min[list], Max[list]}. Over an association, uses\n"
        "\tthe values.");
    symtab_add_builtin("Nearest", builtin_nearest);
    symtab_get_def("Nearest")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Nearest",
        "Nearest[list, x]\n\tGives the element of list closest to x, as a list.\n"
        "\tAll elements tied at the minimum distance Abs[element - x] are\n"
        "\treturned, in their original order; an empty list gives {}.\n"
        "\tReturns unevaluated unless every distance is a real number, so a\n"
        "\tsymbolic element or target leaves the expression unchanged rather\n"
        "\tthan dropping it from the result.");
    symtab_add_builtin("FindClusters", builtin_find_clusters);
    symtab_get_def("FindClusters")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindClusters",
        "FindClusters[list]\n\tPartitions a 1D numeric list into clusters of\n"
        "\tnearby elements, as a list of lists. Clusters appear in order of the\n"
        "\tfirst occurrence of a member; elements keep their input order.\n"
        "FindClusters[list, n]\n\tGives exactly n clusters, capped at the number\n"
        "\tof distinct values.\n"
        "FindClusters[list, UpTo[n]]\n\tGives at most n clusters, and fewer when\n"
        "\tthe data suggests fewer.\n"
        "FindClusters[list, spec, Method -> m]\n\tUses algorithm m: Agglomerate,\n"
        "\tSpanningTree, KMeans, KMedoids, Spectral, DBSCAN, GaussianMixture,\n"
        "\tJarvisPatrick, MeanShift or NeighborhoodContraction. KMeans and\n"
        "\tKMedoids require a count; the density methods require Automatic.\n"
        "\tOptions: Method, DistanceFunction (Automatic, EuclideanDistance,\n"
        "\tManhattanDistance or SquaredEuclideanDistance -- all equivalent in\n"
        "\t1D), CriterionFunction and PerformanceGoal (accepted, no effect).\n"
        "\tReturns unevaluated for a non-numeric element, an empty list, a\n"
        "\tmethod incompatible with the count mode, or a list too large for the\n"
        "\tchosen method (Spectral above 2000 elements, MeanShift and\n"
        "\tNeighborhoodContraction above 4000, both being quadratic).");
    /* Defaults surfaced by Options[FindClusters]. Without this the option
     * machinery reports none, even though four are accepted. */
    {
        Expr* rules[4];
        const char* keys[4];
        keys[0] = SYM_Method; keys[1] = SYM_DistanceFunction;
        keys[2] = SYM_CriterionFunction; keys[3] = SYM_PerformanceGoal;
        for (int i = 0; i < 4; i++) {
            Expr* ra[2] = { expr_new_symbol(keys[i]), expr_new_symbol(SYM_Automatic) };
            rules[i] = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
        }
        symtab_set_options("FindClusters",
                           expr_new_function(expr_new_symbol(SYM_List), rules, 4));
    }
    symtab_add_builtin("ListQ", builtin_listq);
    symtab_add_builtin("VectorQ", builtin_vectorq);
    symtab_add_builtin("MatrixQ", builtin_matrixq);
    symtab_add_builtin("HermitianMatrixQ", builtin_hermitian_matrix_q);
    symtab_add_builtin("SymmetricMatrixQ", builtin_symmetric_matrix_q);
    symtab_add_builtin("SquareMatrixQ", builtin_square_matrix_q);
    symtab_add_builtin("DiagonalMatrixQ", builtin_diagonal_matrix_q);
    symtab_add_builtin("UpperTriangularMatrixQ", builtin_upper_triangular_matrix_q);
    symtab_add_builtin("Subsets", builtin_subsets);
    symtab_get_def("Subsets")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Subsets",
        "Subsets[list]\n\tGives all subsets of list (the power set), ordered by\n"
        "\tincreasing length and lexicographically by element position within\n"
        "\teach length. The head of list is kept on the subsets.\n"
        "Subsets[list, n]\n\tGives subsets of length 0 through n.\n"
        "Subsets[list, {n}]\n\tGives subsets of length exactly n.\n"
        "Subsets[list, {nmin, nmax}]\n\tGives subsets whose length lies in the\n"
        "\tinclusive range nmin to nmax; a third element gives a length step.\n"
        "Subsets[list, spec, s]\n\tGives only the first s subsets spec would\n"
        "\tproduce, generated lazily.");
    symtab_add_builtin("Riffle", builtin_riffle);
    symtab_get_def("Riffle")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Riffle",
        "Riffle[list, x]\n\tInterleaves x into the gaps between successive\n"
        "\telements of list, giving {e1, x, e2, x, ..., x, en}. Nothing is\n"
        "\tplaced before the first or after the last element, so a list of\n"
        "\tlength 0 or 1 comes back unchanged.\n"
        "Riffle[list, {x1, x2, ...}]\n\tUses the xi cyclically, filling the\n"
        "\tn - 1 gaps left to right; separators beyond the last gap are\n"
        "\tunused. The head of list is preserved.");

    symtab_add_builtin("Subdivide", builtin_subdivide);
    symtab_get_def("Subdivide")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Subdivide",
        "Subdivide[n]\n\tGives the list {0, 1/n, 2/n, ..., 1} of n + 1 equally\n"
        "\tspaced points spanning 0 to 1, including both endpoints.\n"
        "Subdivide[max, n]\n\tGives n + 1 equally spaced points spanning 0 to\n"
        "\tmax.\n"
        "Subdivide[min, max, n]\n\tGives n + 1 equally spaced points spanning\n"
        "\tmin to max; point i is min + i (max - min)/n. Descending intervals\n"
        "\t(min > max) are allowed and give a negative step. Exact input gives\n"
        "\texact results in lowest terms, with both endpoints exact. Returns\n"
        "\tunevaluated unless n is a positive integer.");

    symtab_get_def("Table")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_get_def("Range")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Array")->attributes |= ATTR_PROTECTED;
    symtab_get_def("ConstantArray")->attributes |= ATTR_PROTECTED;
    symtab_get_def("ArrayFlatten")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Take")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Drop")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Flatten")->attributes |= ATTR_PROTECTED;
    symtab_get_def("FlattenAt")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Partition")->attributes |= ATTR_PROTECTED;
    symtab_get_def("RotateLeft")->attributes |= ATTR_PROTECTED;
    symtab_get_def("RotateRight")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Reverse")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Rescale")->attributes |= ATTR_NUMERICFUNCTION | ATTR_PROTECTED;
    symtab_get_def("PadRight")->attributes |= ATTR_PROTECTED;
    symtab_get_def("PadLeft")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Transpose")->attributes |= ATTR_PROTECTED;
    symtab_get_def("ConjugateTranspose")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Tally")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Union")->attributes |= ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_get_def("Intersection")->attributes |= ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED;
    symtab_get_def("Complement")->attributes |= ATTR_PROTECTED;
    symtab_get_def("DeleteDuplicates")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Split")->attributes |= ATTR_PROTECTED;
    symtab_get_def("SplitBy")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Total")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Accumulate")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Differences")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Ratios")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Commonest")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Min")->attributes |= ATTR_FLAT | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_ORDERLESS | ATTR_PROTECTED;
    symtab_get_def("Max")->attributes |= ATTR_FLAT | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_ORDERLESS | ATTR_PROTECTED;
    symtab_get_def("ListQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("VectorQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("MatrixQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("HermitianMatrixQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("SymmetricMatrixQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("SquareMatrixQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("DiagonalMatrixQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("UpperTriangularMatrixQ")->attributes |= ATTR_PROTECTED;

    symtab_set_docstring("Total", "Total[list]\n\tgives the total of the elements in list.\nTotal[list, n]\n\ttotals all elements down to level n.\nTotal[list, {n}]\n\ttotals elements at level n.\nTotal[list, {n1, n2}]\n\ttotals elements at levels n1 through n2.");
}
