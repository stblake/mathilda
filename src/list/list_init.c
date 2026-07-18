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
    symtab_add_builtin("Partition", builtin_partition);
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
    symtab_add_builtin("DeleteDuplicates", builtin_deleteduplicates);
    symtab_add_builtin("DeleteDuplicatesBy", builtin_deleteduplicatesby);
    symtab_get_def("DeleteDuplicatesBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DeleteDuplicatesBy",
        "DeleteDuplicatesBy[expr, f]\n\tKeeps the first element for each distinct\n"
        "\tf[element], preserving order. Over an association, f is applied to the\n"
        "\tvalues and the surviving entries are kept (keys preserved).");
    symtab_add_builtin("Split", builtin_split);
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
    symtab_add_builtin("ListQ", builtin_listq);
    symtab_add_builtin("VectorQ", builtin_vectorq);
    symtab_add_builtin("MatrixQ", builtin_matrixq);
    symtab_add_builtin("HermitianMatrixQ", builtin_hermitian_matrix_q);
    symtab_add_builtin("SymmetricMatrixQ", builtin_symmetric_matrix_q);
    symtab_add_builtin("SquareMatrixQ", builtin_square_matrix_q);
    symtab_add_builtin("DiagonalMatrixQ", builtin_diagonal_matrix_q);
    symtab_add_builtin("UpperTriangularMatrixQ", builtin_upper_triangular_matrix_q);

    symtab_get_def("Table")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_get_def("Range")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Array")->attributes |= ATTR_PROTECTED;
    symtab_get_def("ConstantArray")->attributes |= ATTR_PROTECTED;
    symtab_get_def("ArrayFlatten")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Take")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Drop")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Flatten")->attributes |= ATTR_PROTECTED;
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
    symtab_get_def("DeleteDuplicates")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Split")->attributes |= ATTR_PROTECTED;
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
