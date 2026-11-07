/* linalg_init - the linear-algebra module entry point.
 *
 * Registers the per-builtin symbols, attributes, and docstrings.  Each
 * builtin lives in its own translation unit in src/linalg/ (dot.c, det.c,
 * cross.c, ...).
 *
 * The neighbouring modules Inverse / PseudoInverse (matinv.c),
 * RowReduce / LinearSolve (matsol.c), least squares (matlstsq.c), and
 * Eigenvalues / Eigenvectors (mateigen.c) register themselves through
 * their own init functions, called from core_init().
 */

#include "linalg.h"
#include "hnf.h"
#include "symtab.h"
#include "attr.h"

void linalg_init(void) {
    symtab_add_builtin("Dot", builtin_dot);
    symtab_get_def("Dot")->attributes |= ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED;
    symtab_add_builtin("Det", builtin_det);
    symtab_get_def("Det")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Cross", builtin_cross);
    symtab_get_def("Cross")->attributes |= ATTR_PROTECTED;
    rotationmatrix_init();
    symtab_add_builtin("Projection", builtin_projection);
    symtab_get_def("Projection")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("TensorProduct", builtin_tensorproduct);
    symtab_get_def("TensorProduct")->attributes |= ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED;
    symtab_add_builtin("Norm", builtin_norm);
    symtab_get_def("Norm")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Normalize", builtin_normalize);
    symtab_get_def("Normalize")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Tr", builtin_tr);
    symtab_get_def("Tr")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Diagonal", builtin_diagonal);
    symtab_get_def("Diagonal")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("IdentityMatrix", builtin_identitymatrix);
    symtab_get_def("IdentityMatrix")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("DiagonalMatrix", builtin_diagonalmatrix);
    symtab_get_def("DiagonalMatrix")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("KroneckerProduct", builtin_kroneckerproduct);
    symtab_get_def("KroneckerProduct")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KroneckerProduct",
        "KroneckerProduct[a, b, ...] gives the Kronecker product of the arrays a, b, ...\n"
        "For matrices it is the block matrix of a[[i, j]] b; for two vectors the outer\n"
        "product (a column times a row), and a vector after a matrix is a row. Arrays of\n"
        "any rank are accepted; with three or more arguments all but the last must have\n"
        "rank 1 or 2. Machine operands give a packed result; integer products stay exact.");
    symtab_add_builtin("UpperTriangularize", builtin_uppertriangularize);
    symtab_get_def("UpperTriangularize")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("UpperTriangularize",
        "UpperTriangularize[m] gives m with every entry below the leading diagonal replaced by 0.\n"
        "UpperTriangularize[m, k] keeps the entries on and above the k-th diagonal (k > 0 above,\n"
        "k < 0 below the leading one). m may be rectangular; a matrix holding machine numbers\n"
        "gets the zero 0.");
    symtab_add_builtin("LowerTriangularize", builtin_lowertriangularize);
    symtab_get_def("LowerTriangularize")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LowerTriangularize",
        "LowerTriangularize[m] gives m with every entry above the leading diagonal replaced by 0.\n"
        "LowerTriangularize[m, k] keeps the entries on and below the k-th diagonal (k > 0 above,\n"
        "k < 0 below the leading one). m may be rectangular; a matrix holding machine numbers\n"
        "gets the zero 0.");
    symtab_add_builtin("HilbertMatrix", builtin_hilbertmatrix);
    symtab_get_def("HilbertMatrix")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("HankelMatrix", builtin_hankelmatrix);
    symtab_get_def("HankelMatrix")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ToeplitzMatrix", builtin_toeplitzmatrix);
    symtab_get_def("ToeplitzMatrix")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("VandermondeMatrix", builtin_vandermondematrix);
    symtab_get_def("VandermondeMatrix")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("MatrixPower", builtin_matrixpower);
    symtab_get_def("MatrixPower")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("PositiveDefiniteMatrixQ",
                       builtin_positive_definite_matrix_q);
    symtab_get_def("PositiveDefiniteMatrixQ")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("NegativeDefiniteMatrixQ",
                       builtin_negative_definite_matrix_q);
    symtab_get_def("NegativeDefiniteMatrixQ")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("LatticeReduce", builtin_latticereduce);
    symtab_get_def("LatticeReduce")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("FindIntegerNullVector", builtin_findintegernullvector);
    symtab_get_def("FindIntegerNullVector")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("HermiteDecomposition", builtin_hermite_decomposition);
    symtab_get_def("HermiteDecomposition")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("HermiteDecomposition",
        "HermiteDecomposition[m]\n"
        "\tGives the Hermite normal form decomposition {u, r} of the integer\n"
        "\tmatrix m: u is unimodular (Abs[Det[u]] == 1), r is the row Hermite\n"
        "\tnormal form, and u . m == r.  r is in echelon shape with positive\n"
        "\tpivots and entries above each pivot reduced into [0, pivot).");
}
