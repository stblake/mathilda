/*
 * sym_intern.h
 *
 * Global string interner for symbol names. Returns a stable, canonical
 * pointer for each unique name string -- two symbol expressions with
 * the same name share the same `const char*`, so identity tests reduce
 * to pointer equality.
 *
 * Lifetime: an interned pointer is valid until intern_clear() is called
 * (typically at program shutdown). The interner owns the string
 * memory; callers must NOT free returned pointers.
 *
 * The interner lazily initializes on the first call to intern_symbol(),
 * so it is safe to call before any explicit init step.
 */

#ifndef SYM_INTERN_H
#define SYM_INTERN_H

#ifdef __cplusplus
extern "C" {
#endif

/* Return a canonical, table-owned pointer for `name`. Equal name strings
 * always return the same pointer. NULL input returns NULL. */
const char* intern_symbol(const char* name);

/* Flag every name interned so far as a System symbol. Called once, after all
 * System symbols have been interned (end of the kernel's C-side init) and
 * before any user/package source is parsed. Names interned afterwards (a
 * package's private symbols, user identifiers) are NOT flagged. */
void intern_mark_all_system(void);

/* True if `name` was flagged by intern_mark_all_system() -- i.e. it is a
 * built-in System symbol. The context resolver uses this to resolve such a
 * bare name to System` instead of qualifying it into the current context,
 * even when the symbol has no SymbolDef yet (pure structural heads and
 * constants like List, Rule, Automatic, Infinity are created lazily). */
int intern_is_system(const char* name);

/* Free every interned string and clear the table. Any previously-returned
 * pointer becomes invalid. Intended for end-of-program cleanup. */
void intern_clear(void);

#ifdef __cplusplus
}
#endif

#endif /* SYM_INTERN_H */
