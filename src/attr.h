#ifndef ATTR_H
#define ATTR_H

#include "expr.h"
#include <stdint.h>

// Attribute bitmasks
#define ATTR_NONE           0
#define ATTR_HOLDFIRST      (1 << 0)
#define ATTR_HOLDREST       (1 << 1)
#define ATTR_HOLDALL        (ATTR_HOLDFIRST | ATTR_HOLDREST)
#define ATTR_HOLDALLCOMPLETE (1 << 2)
#define ATTR_FLAT           (1 << 3)
#define ATTR_ORDERLESS      (1 << 4)
#define ATTR_LISTABLE       (1 << 5)
#define ATTR_NUMERICFUNCTION (1 << 6)
#define ATTR_PROTECTED      (1 << 7)
#define ATTR_ONEIDENTITY    (1 << 8)
#define ATTR_NHOLDREST      (1 << 9)
#define ATTR_NHOLDFIRST     (1 << 15)
#define ATTR_NHOLDALL       (ATTR_NHOLDFIRST | ATTR_NHOLDREST)
#define ATTR_LOCKED          (1 << 10)
#define ATTR_READPROTECTED   (1 << 11)
#define ATTR_TEMPORARY       (1 << 12)
#define ATTR_SEQUENCEHOLD    (1 << 13)
#define ATTR_CONSTANT        (1 << 14)

// Get attributes for a given symbol
uint32_t get_attributes(const char* symbol_name);

// Phase 3a: def-threaded fast path -- effective attributes straight from an
// already-resolved SymbolDef, with no symbol-table lookup. Returns ATTR_NONE
// for a NULL def.
struct SymbolDef;
uint32_t get_attributes_def(struct SymbolDef* def);

// Set attributes for a given symbol
void set_attributes(const char* symbol_name, uint32_t attrs);

// Built-in functions
Expr* builtin_attributes(Expr* res);
Expr* builtin_set_attributes(Expr* res);
Expr* builtin_clear_attributes(Expr* res);

// Initialize attributes
void attr_init(void);

#endif // ATTR_H
