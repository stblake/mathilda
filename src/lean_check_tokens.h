/*
 * lean_check_tokens.h -- the command-leading words of Lean 4's parser, for
 * LeanVerify's source check (src/lean_check.c).
 *
 * DERIVED FROM LEAN, NOT WRITTEN BY HAND.  tools/lean_command_tokens.lean
 * reads the `command` parser category from Lean.Parser.parserExtension and
 * prints the first tokens of every leading command parser.  The words below
 * are the [A-Za-z0-9_] runs of every such token that is not already refused
 * by the lexical surface rules (a token containing #, %, $, @[, /- , ... can
 * never reach the parser through an accepted source).  LEAN_COMMAND_WORDS_CORE
 * comes from `import Lean` (a superset of the core Init a LeanVerify source
 * sees), LEAN_COMMAND_WORDS_MATHLIB from `import Mathlib` minus the core list.
 *
 * Kept in sync by tests (ctest):
 *   lean_check_tokens_sync     re-runs the script against the installed core
 *                              Lean and fails if any command-leading token
 *                              is not refused by Lean`Internal`CheckSource;
 *   lean_check_tokens_mathlib  the same with `import Mathlib`, in the Lake
 *                              project of MATHILDA_LEAN_PROJECT.
 * Regenerate with:  sh tests/scripts/lean_check_tokens.sh --words ./Mathilda
 * (Lean v4.35.0-rc3; Mathlib b63f6e8a68d2).
 */
#ifndef MATHILDA_LEAN_CHECK_TOKENS_H
#define MATHILDA_LEAN_CHECK_TOKENS_H

static const char* const LEAN_COMMAND_WORDS_CORE[] = {
    "_tactic", "abbrev", "assert_not_exists", "assert_not_imported", "attribute", "axiom",
    "binder_predicate", "builtin_cbv_simproc", "builtin_cbv_simproc_decl",
    "builtin_dsimproc", "builtin_dsimproc_decl", "builtin_grind_propagator",
    "builtin_initialize", "builtin_simproc", "builtin_simproc_decl", "cbv_simproc",
    "cbv_simproc_decl", "class", "coinductive", "declare_bitwise_int_theorems",
    "declare_bitwise_uint_theorems", "declare_command_config_elab",
    "declare_command_config_elab_legacy", "declare_config_elab",
    "declare_config_elab_legacy", "declare_core_config_elab", "declare_eval_bin",
    "declare_eval_bin_bitwise", "declare_eval_bin_bool_pred", "declare_int_theorems",
    "declare_simp_like_tactic", "declare_sint_simprocs", "declare_syntax_cat",
    "declare_term_config_elab", "declare_uint_simprocs", "declare_uint_theorems", "def",
    "def_eval_config_item", "deprecated_module", "deprecated_syntax", "deriving",
    "docs_to_verso", "dsimproc", "dsimproc_decl", "elab", "elab_rules", "elab_stx_quot",
    "end", "example", "export", "grind_annotated", "grind_pattern", "grind_propagator",
    "import", "include", "inductive", "infix", "infixl", "infixr", "init_grind_norm",
    "init_quot", "initialize", "instance", "local", "macro", "macro_rules", "meta",
    "mutual", "namespace", "noncomputable", "nonrec", "norm_cast_add_elim", "notation",
    "omit", "opaque", "open", "partial", "postfix", "prefix", "private", "protected",
    "public", "recall", "recommended_spelling", "register_builtin_option",
    "register_error_explanation", "register_grind_attr", "register_label_attr",
    "register_linter_set", "register_option", "register_simp_attr", "register_sym_dsimp",
    "register_sym_simp", "register_sym_simp_attr", "register_tactic_tag", "register_try",
    "reprove", "run_cmd", "run_elab", "run_meta", "scoped", "seal", "section",
    "set_library_suggestions", "set_option", "show_panel_widgets", "simproc",
    "simproc_decl", "structure", "syntax", "tactic_extension", "test_extern", "theorem",
    "unif_hint", "universe", "unlock_limits", "unsafe", "unseal", "variable",
    "with_weak_namespace",
    NULL
};

static const char* const LEAN_COMMAND_WORDS_MATHLIB[] = {
    "add_aesop_rules", "alias", "assert_no_sorry", "declare_aesop_rule_sets", "def_wanted",
    "deprecate", "erase_aesop_rules", "extend_docs", "guard_min_heartbeats",
    "initialize_simps_projections", "insert_to_additive_translation",
    "insert_to_dual_translation", "instance_wanted", "irreducible_def", "lemma",
    "library_note", "lrat_proof", "mk_iff_of_inductive_prop", "name_poly_vars",
    "name_power_vars", "notation3", "proof_wanted", "register_hint", "sudo",
    "suppress_compilation", "theorem_wanted", "to_additive_name_hint",
    "to_dual_insert_cast", "to_dual_insert_cast_fun", "to_dual_name_hint", "unset_option",
    "unsuppress_compilation", "variables", "whatsnew",
    NULL
};

/* Tokens of Lean's token table that END in a word character [A-Za-z0-9_] but
 * start (or continue) with other characters, so that Lean's next token can
 * begin inside one of the checker's words: in Mathlib, U+2192 followed by
 * `otheorem` is the token U+2192 `o` (the order-homomorphism arrow) followed
 * by the keyword `theorem`, where the checker sees the word `otheorem`.  Wherever one of these
 * occurs, the checker also checks the word that follows it.  Only tokens
 * inside the lexical surface are listed (the others never reach the word
 * check); core Lean has none.  From `import Mathlib`; kept in sync by
 * lean_check_tokens_mathlib. */
static const char* const LEAN_GLUED_TOKENS_MATHLIB[] = {
    "<i", "<m", "[MOD", "[PMOD", "[SMOD", "[ZMOD", "~r", "\xCE\xB2_", "\xCE\xB5_",
    "\xCE\xB7_", "\xE2\x84\x9A\xE2\x89\xA5" "0", "\xE2\x86\x92*o", "\xE2\x86\x92+*o",
    "\xE2\x86\x92+o", "\xE2\x86\x92" "CO", "\xE2\x86\x92" "Co", "\xE2\x86\x92g",
    "\xE2\x86\x92o", "\xE2\x86\x92r", "\xE2\x88\x98L", "\xE2\x88\x98SL", "\xE2\x88\xA3_",
    "\xE2\x89\xA4i", "\xE2\x89\xA4m", "\xE2\x89\xA4s",
    NULL
};

#endif /* MATHILDA_LEAN_CHECK_TOKENS_H */
