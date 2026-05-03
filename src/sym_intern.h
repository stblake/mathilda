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

/* Free every interned string and clear the table. Any previously-returned
 * pointer becomes invalid. Intended for end-of-program cleanup. */
void intern_clear(void);

#ifdef __cplusplus
}
#endif

#endif /* SYM_INTERN_H */
