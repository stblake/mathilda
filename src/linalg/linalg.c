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
#include "symtab.h"
#include "attr.h"

void linalg_init(void) {
    symtab_add_builtin("Dot", builtin_dot);
    symtab_get_def("Dot")->attributes |= ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED;
    symtab_add_builtin("Det", builtin_det);
    symtab_get_def("Det")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Cross", builtin_cross);
    symtab_get_def("Cross")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Norm", builtin_norm);
    symtab_get_def("Norm")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Normalize", builtin_normalize);
    symtab_get_def("Normalize")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Tr", builtin_tr);
    symtab_get_def("Tr")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("IdentityMatrix", builtin_identitymatrix);
    symtab_get_def("IdentityMatrix")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("DiagonalMatrix", builtin_diagonalmatrix);
    symtab_get_def("DiagonalMatrix")->attributes |= ATTR_PROTECTED;
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
}
