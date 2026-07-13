/*
 * sym_intern.c
 *
 * Phase 2 (EVAL_SYMTAB_IMPROVEMENTS): the string interner has been MERGED into
 * the symbol table. The unified hash table in symtab.c is now keyed on the
 * symbol name and backs both interning (canonical `const char*`) and definition
 * lookup, so resolving a symbol no longer hashes twice across two separate
 * tables.
 *
 * The sym_intern.h API (intern_symbol, intern_is_system, intern_mark_all_system,
 * intern_clear) is therefore implemented in symtab.c over that unified table.
 * This translation unit is intentionally empty apart from this note; the
 * typedef below keeps it from being an empty ISO C translation unit.
 */

typedef int mathilda_sym_intern_merged_into_symtab;
