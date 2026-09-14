/*
 * readlist.h — ReadList[] file-reading builtin (src/io subsystem).
 *
 * ReadList reads the objects contained in a file and returns them as a
 * List, directed by a read type (Expression, Number, Real, Word, Record,
 * String, Character, Byte) and modified by separator options.  This is the
 * first "read objects from a file" builtin in the tree; there is no stream
 * layer, so ReadList always opens the named file, reads it, and closes it.
 */
#ifndef MATHILDA_IO_READLIST_H
#define MATHILDA_IO_READLIST_H

#include "expr.h"

/* ReadList["file"], ReadList["file", type], ReadList["file", {t1,...}],
 * ReadList["file", types, n], each with trailing separator options. */
Expr* builtin_readlist(Expr* res);

/* Registers ReadList, its Protected attribute, and its default Options. */
void readlist_init(void);

#endif /* MATHILDA_IO_READLIST_H */
