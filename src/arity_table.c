/* Mathilda — the declared argument counts of the builtins (see arity.h).
 *
 * One row per head: the reference's positional argument count [min, max]
 * (max -1: unbounded), and, for a head that takes options, the first index
 * at which option names are checked against Options[head] (OPT_NOCHECK: the
 * head takes options, but Options[head] is not the full list it reads).  A
 * head with a non-empty Options[head] takes options whatever its row says.
 *
 * A row is a claim about every call of the head: a max one too small turns a
 * valid call into a message.  So a row is added only with the head's
 * documented signatures in hand (docs/spec/builtins/), operator forms
 * included -- Map[f] is a call with one argument.
 *
 * Applied by arity_table_init, at the end of core_init, once every subsystem
 * has registered its builtins and options.
 */

#include "arity.h"

#define ANY          (-1)    /* max: unbounded */
#define NOOPT        (-2)    /* takes no options */
#define OPT_NOCHECK  (-1)    /* takes options; names not checked */

typedef struct {
    const char* name;
    short min, max;
    short vfrom;    /* NOOPT, OPT_NOCHECK, or the first index whose names are checked */
} ArityRow;

static const ArityRow ARITY_ROWS[] = {
    /* ---- evaluation control, scoping, assignment ---- */
    { "If",               2, 4, NOOPT },
    { "For",              3, 4, NOOPT },
    { "While",            1, 2, NOOPT },
    { "Catch",            1, 3, NOOPT },
    { "Module",           2, 2, NOOPT },
    { "Goto",             1, 1, NOOPT },
    { "AbortProtect",     1, 1, NOOPT },
    { "CheckAbort",       2, 2, NOOPT },
    { "TimeConstrained",  2, 3, NOOPT },
    { "Quiet",            1, 3, NOOPT },
    { "Timing",           1, 1, NOOPT },
    { "AbsoluteTiming",   1, 1, NOOPT },
    { "RepeatedTiming",   1, 2, NOOPT },
    { "Pause",            1, 1, NOOPT },
    { "Monitor",          2, 2, NOOPT },
    { "Assuming",         2, 2, NOOPT },
    { "Dialog",           0, 1, OPT_NOCHECK },
    { "SessionTime",      0, 0, NOOPT },
    { "TimeUsed",         0, 0, NOOPT },
    { "Set",              2, 2, NOOPT },
    { "SetDelayed",       2, 2, NOOPT },
    { "Rule",             2, 2, NOOPT },
    { "RuleDelayed",      2, 2, NOOPT },
    { "AddTo",            2, 2, NOOPT },
    { "SubtractFrom",     2, 2, NOOPT },
    { "TimesBy",          2, 2, NOOPT },
    { "DivideBy",         2, 2, NOOPT },
    { "AppendTo",         2, 2, NOOPT },
    { "PrependTo",        2, 2, NOOPT },
    { "ApplyTo",          2, 2, NOOPT },
    { "Increment",        1, 1, NOOPT },
    { "Decrement",        1, 1, NOOPT },
    { "PreIncrement",     1, 1, NOOPT },
    { "PreDecrement",     1, 1, NOOPT },
    { "SetAttributes",    2, 2, NOOPT },
    { "ClearAttributes",  2, 2, NOOPT },
    { "DownValues",       1, 1, OPT_NOCHECK },
    { "UpValues",         1, 1, OPT_NOCHECK },
    { "Activate",         1, 2, NOOPT },
    { "Block",            2, 2, NOOPT },
    { "Unset",            1, 1, NOOPT },
    { "Attributes",       1, 1, NOOPT },
    { "Options",          1, 2, NOOPT },
    { "BeginPackage",     1, 2, NOOPT },
    { "MaxMemoryUsed",    0, 1, NOOPT },

    /* ---- arithmetic and numeric functions ---- */
    { "Minus",            1, 1, NOOPT },
    { "Subtract",         2, 2, NOOPT },
    { "Divide",           2, 2, NOOPT },
    { "Abs",              1, 1, NOOPT },
    { "Arg",              1, 1, NOOPT },
    { "Sign",             1, 1, NOOPT },
    { "ReIm",             1, 1, NOOPT },
    { "Ramp",             1, 1, NOOPT },
    { "Boole",            1, 1, NOOPT },
    { "IntegerPart",      1, 1, NOOPT },
    { "Round",            1, 2, NOOPT },
    { "Rationalize",      1, 2, NOOPT },
    { "Clip",             1, 3, NOOPT },
    { "Factorial",        1, 1, NOOPT },
    { "FactorialPower",   2, 3, NOOPT },
    { "Gamma",            1, 3, NOOPT },
    { "BesselJ",          2, 2, NOOPT },
    { "BesselJZero",      2, 3, NOOPT },
    { "AiryBiPrimeZero",  1, 2, NOOPT },
    { "AngerJ",           2, 3, NOOPT },
    { "WeberE",           2, 3, NOOPT },
    { "KelvinBer",        1, 2, NOOPT },
    { "KelvinBei",        1, 2, NOOPT },
    { "KelvinKer",        1, 2, NOOPT },
    { "KelvinKei",        1, 2, NOOPT },
    { "StieltjesGamma",   1, 2, NOOPT },
    { "Piecewise",        1, 2, NOOPT },
    { "MathieuCharacteristicA", 2, 2, NOOPT },
    { "QPochhammer",      1, 3, NOOPT },
    { "Precision",        1, 1, NOOPT },
    { "Accuracy",         1, 1, NOOPT },
    { "SetAccuracy",      2, 2, NOOPT },
    { "FromContinuedFraction", 1, 1, NOOPT },
    { "Signature",        1, 1, NOOPT },
    { "PermutationCycles", 1, 2, NOOPT },
    { "GroupElements",    1, 2, NOOPT },

    /* ---- predicates ---- */
    { "AtomQ",            1, 1, NOOPT },
    { "BooleanQ",         1, 1, NOOPT },
    { "EvenQ",            1, 1, NOOPT },
    { "OddQ",             1, 1, NOOPT },
    { "IntegerQ",         1, 1, NOOPT },
    { "NumberQ",          1, 1, NOOPT },
    { "NumericQ",         1, 1, NOOPT },
    { "ExactNumberQ",     1, 1, NOOPT },
    { "InexactNumberQ",   1, 1, NOOPT },
    { "MachineNumberQ",   1, 1, NOOPT },
    { "StringQ",          1, 1, NOOPT },
    { "LetterQ",          1, 1, NOOPT },
    { "ListQ",            1, 1, NOOPT },
    { "MissingQ",         1, 1, NOOPT },
    { "TrueQ",            1, 1, NOOPT },
    { "OptionQ",          1, 1, NOOPT },
    { "Not",              1, 1, NOOPT },
    { "Implies",          2, 2, NOOPT },
    { "MatchQ",           1, 2, NOOPT },
    { "MatrixQ",          1, 2, NOOPT },
    { "OrderedQ",         1, 2, NOOPT },
    { "IntervalMemberQ",  2, 2, NOOPT },
    { "VectorQ",          1, 2, NOOPT },
    { "AllTrue",          1, 3, NOOPT },   /* AllTrue[test] operator form */
    { "AnyTrue",          1, 3, NOOPT },
    { "NoneTrue",         1, 3, NOOPT },
    { "LeafCount",        1, 1, OPT_NOCHECK },   /* Heads -> True */

    /* ---- structural and list functions ---- */
    { "Keys",             1, 2, NOOPT },
    { "Values",           1, 2, NOOPT },
    { "Length",           1, 1, NOOPT },
    { "Head",             1, 2, NOOPT },
    { "First",            1, 2, NOOPT },
    { "Last",             1, 2, NOOPT },
    { "Most",             1, 1, NOOPT },
    { "Rest",             1, 1, NOOPT },
    { "Append",           1, 2, NOOPT },
    { "Prepend",          1, 2, NOOPT },
    { "Insert",           2, 3, NOOPT },
    { "Delete",           1, 2, NOOPT },
    { "Reverse",          1, 2, NOOPT },
    { "RotateLeft",       1, 2, NOOPT },
    { "RotateRight",      1, 2, NOOPT },
    { "Sort",             1, 2, NOOPT },
    { "Split",            1, 2, NOOPT },
    { "SplitBy",          2, 2, NOOPT },
    { "Partition",        2, 6, NOOPT },
    { "Tally",            1, 2, NOOPT },
    { "DeleteDuplicates", 1, 2, NOOPT },
    { "CountDistinct",    1, 2, NOOPT },
    { "Accumulate",       1, 1, OPT_NOCHECK },   /* Method -> "CompensatedSummation" */
    { "TakeWhile",        2, 2, NOOPT },
    { "LengthWhile",      2, 2, NOOPT },
    { "FlattenAt",        2, 2, NOOPT },
    { "MapAt",            2, 3, NOOPT },
    { "Thread",           1, 3, NOOPT },
    { "Distribute",       1, 5, NOOPT },
    { "NestList",         3, 3, NOOPT },
    { "Nest",             3, 3, NOOPT },
    { "Fold",             1, 3, NOOPT },   /* Fold[f] operator form */
    { "FoldList",         1, 3, NOOPT },
    { "Position",         1, 4, OPT_NOCHECK },   /* Position[pattern] operator form */
    { "ReverseSort",      1, 2, NOOPT },
    { "ToString",         1, 2, OPT_NOCHECK },
    { "NestWhile",        3, 6, NOOPT },
    { "NestWhileList",    3, 6, NOOPT },
    { "FirstPosition",    2, 4, OPT_NOCHECK },
    { "Map",              1, 3, 2 },
    { "Count",            1, 3, 2 },    /* Count[pattern] operator form */
    { "MemberQ",          1, 3, 2 },
    { "FreeQ",            1, 3, 2 },
    { "Level",            2, 3, 3 },    /* Level[expr, spec, h] */
    { "Total",            1, 2, OPT_NOCHECK },
    { "MapAll",           1, 2, 2 },
    { "Apply",            1, 3, 2 },
    { "Scan",             1, 3, 2 },
    { "Replace",          1, 3, 2 },
    { "ReplaceAll",       1, 2, NOOPT },
    { "ReplaceList",      2, 3, NOOPT },
    { "Dispatch",         1, 1, NOOPT },
    { "ByteCount",        1, 1, NOOPT },
    { "Hash",             1, 3, NOOPT },
    { "FromPackedArray",  1, 1, NOOPT },
    { "PositionIndex",    1, 1, NOOPT },
    { "KeyComplement",    1, 1, NOOPT },
    { "KeyValueMap",      1, 2, NOOPT },
    { "Merge",            1, 2, NOOPT },
    { "GroupBy",          1, 3, NOOPT },
    { "CountsBy",         1, 2, NOOPT },
    { "IntersectingQ",    2, 2, OPT_NOCHECK },

    /* ---- strings, symbols, files ---- */
    { "StringLength",     1, 1, NOOPT },
    { "Decapitalize",     1, 1, NOOPT },
    { "SymbolName",       1, 1, NOOPT },
    { "Environment",      1, 1, NOOPT },
    { "OpenRead",         1, 1, OPT_NOCHECK },
    { "FileNameSplit",    1, 1, OPT_NOCHECK },
    { "StringContainsQ",  1, 2, 2 },
    { "StringCount",      2, 2, 2 },
    { "StringDelete",     1, 2, 2 },
    { "StringMatchQ",     1, 2, 2 },

    /* ---- statistics ---- */
    { "Mean",             1, 1, NOOPT },
    { "Covariance",       1, 2, NOOPT },
    { "CentralMoment",    2, 2, NOOPT },

    /* ---- algebra ---- */
    { "FunctionExpand",   1, 2, OPT_NOCHECK },
    { "TrigToExp",        1, 1, NOOPT },
    { "ExpToTrig",        1, 1, NOOPT },
    { "Variables",        1, 1, NOOPT },
    { "Refine",           1, 2, OPT_NOCHECK },
    { "TrigReduce",       1, 1, NOOPT },
    { "ExpandDenominator", 1, 1, NOOPT },
    { "Collect",          2, 3, OPT_NOCHECK },
    { "FactorTermsList",  1, 2, NOOPT },
    { "FromCoefficientRules", 2, 2, NOOPT },
    { "Subresultants",    3, 3, OPT_NOCHECK },
    { "SubresultantPolynomials", 3, 3, OPT_NOCHECK },
    { "SubresultantPolynomialRemainders", 3, 3, OPT_NOCHECK },
    { "AlgebraicNumberNorm", 1, 1, OPT_NOCHECK },
    { "ToNumberField",    1, 2, NOOPT },
    { "NumberFieldIntegralBasis", 1, 1, NOOPT },
    { "RootReduce",       1, 1, 1 },
    { "FactorList",       1, 1, 1 },
    { "Cancel",           1, 1, 1 },
    { "PolynomialGCD",    1, ANY, 1 },
    { "PolynomialLCM",    1, ANY, 1 },
    { "InterpolatingPolynomial", 2, 2, 2 },
    { "PolynomialReduce", 2, 3, 2 },
    { "GroebnerBasis",    1, 3, 2 },

    /* ---- linear algebra ---- */
    { "Det",              1, 1, OPT_NOCHECK },
    { "Norm",             1, 2, NOOPT },
    { "Tr",               1, 3, NOOPT },
    { "Eigensystem",      1, 2, 1 },
    { "Eigenvalues",      1, 2, 1 },
    { "Eigenvectors",     1, 2, 1 },
    { "LUDecomposition",  1, 1, NOOPT },
    { "QRDecomposition",  1, 1, 1 },
    { "SchurDecomposition", 1, 1, 1 },
    { "MatrixFunction",   2, 2, 2 },
    { "LinearSolve",      1, 2, 2 },
    { "VandermondeMatrix", 1, 2, NOOPT },
    { "HilbertMatrix",    1, 1, 1 },
    { "UnitVector",       1, 2, 1 },
    { "DiagonalMatrixQ",  1, 2, 1 },
    { "UpperTriangularMatrixQ", 1, 2, 1 },
    { "LowerTriangularMatrixQ", 1, 2, 1 },
    { "Inverse",          1, 1, OPT_NOCHECK },
    { "NullSpace",        1, 1, OPT_NOCHECK },
    { "IdentityMatrix",   1, 2, OPT_NOCHECK },
    { "SymmetricMatrixQ", 1, 1, OPT_NOCHECK },
    { "Inner",            3, 5, NOOPT },

    /* ---- graphs ---- */
    { "GraphPlot",        1, 1, OPT_NOCHECK },

    /* ---- calculus and numerics ---- */
    { "Limit",            2, 2, 2 },
    { "FindMinimum",      2, ANY, 2 },
    { "FindMaximum",      2, ANY, 2 },
    { "NMinimize",        2, 3, 2 },
    { "NMaximize",        2, 3, 2 },
    { "NSum",             2, ANY, 2 },
    { "NLimit",           2, 2, 2 },
    { "NSeries",          2, 2, 2 },
    { "NDSolve",          2, ANY, 2 },
    { "NDSolveValue",     2, ANY, 2 },
    { "DSolve",           2, 3, 3 },
    { "Maximize",         2, 3, 2 },
    { "Minimize",         2, 3, 2 },
    { "Fit",              1, 4, 1 },   /* Fit[{m, v}] is the design-matrix form */

    /* ---- images ---- */
    { "AlphaChannel",     1, 1, NOOPT },
    { "Erosion",          2, 2, OPT_NOCHECK },
    { "FindThreshold",    1, 1, OPT_NOCHECK },
    { "LocalAdaptiveBinarize", 2, 3, OPT_NOCHECK },
    { "RandomImage",      0, 2, 0 },
    { "WaveletPhi",       1, 2, 2 },

    /* ---- graphics (Wave 5 G3): a non-option argument past the iterators
     * is General::nonopt, not a silent unevaluated call (audit Plot#6,
     * ParametricPlot#6, ListPlot#3).  Option names are not checked here:
     * the plotters pass Graphics options through. ---- */
    { "Plot",             2, 2, OPT_NOCHECK },
    { "ParametricPlot",   2, 3, OPT_NOCHECK },
    { "PolarPlot",        2, 2, OPT_NOCHECK },
    { "Plot3D",           3, 3, OPT_NOCHECK },
    { "ParametricPlot3D", 2, 3, OPT_NOCHECK },
    { "ListPlot",         1, 1, OPT_NOCHECK },
    /* ---- Wave 5 G2 (structural, data structures, lists) ---- */
    { "Expand",           1, 2, OPT_NOCHECK },   /* Expand[expr, patt]; Modulus, Trig */
    { "Factor",           1, 1, OPT_NOCHECK },
    { "FactorTerms",      1, 2, OPT_NOCHECK },   /* Modulus */
    { "PolynomialQ",      1, 2, NOOPT },
    { "Permutations",     1, 2, NOOPT },
    { "Discard",          1, 3, NOOPT },         /* Discard[crit] operator form */
    { "AssociationThread", 1, 2, NOOPT },
    { "AssociationComap", 1, 2, NOOPT },

    /* ---- files and packed arrays ---- */
    { "ReadList",         1, 3, 1 },
    { "ToNDArray",        1, 1, 1 },

    /* ---- Wave 5 G7 (arithmetic, strings, random): argument-count messages ---- */
    { "N",                1, 2, NOOPT },
    { "Re",               1, 1, NOOPT },
    { "Im",               1, 1, NOOPT },
    { "Quotient",         2, 3, NOOPT },
    { "LucasL",           1, 2, NOOPT },
    { "Capitalize",       1, 2, NOOPT },   /* Capitalize[s, "TitleCase"] */
    { "CharacterRange",   2, 2, NOOPT },
    { "DigitQ",           1, 1, NOOPT },
    { "LowerCaseQ",       1, 1, NOOPT },
    { "UpperCaseQ",       1, 1, NOOPT },
    { "ToLowerCase",      1, 1, NOOPT },
    { "ToUpperCase",      1, 1, NOOPT },
    { "ToCharacterCode",  1, 2, NOOPT },
    { "RandomChoice",     1, 2, NOOPT },
    { "RandomInteger",    0, 2, NOOPT },
};

/* Heads whose choice of two counts the reference words as a range:
 * "While called with 3 arguments; between 1 and 2 arguments are expected." */
static const char* const BETWEEN_ROWS[] = { "For", "While", "IdentityMatrix" };

void arity_table_init(void) {
    for (size_t i = 0; i < sizeof(ARITY_ROWS) / sizeof(ARITY_ROWS[0]); i++) {
        const ArityRow* r = &ARITY_ROWS[i];
        if (r->vfrom == NOOPT) symtab_set_arity(r->name, r->min, r->max);
        else symtab_set_arity_opts(r->name, r->min, r->max, r->vfrom);
    }
    for (size_t i = 0; i < sizeof(BETWEEN_ROWS) / sizeof(BETWEEN_ROWS[0]); i++)
        symtab_set_arity_between(BETWEEN_ROWS[i]);
}
