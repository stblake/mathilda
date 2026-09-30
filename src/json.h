#ifndef MATHILDA_JSON_H
#define MATHILDA_JSON_H

/* ---------------------------------------------------------------------------
 * JSON reader and writer (src/json.c) behind ImportString / ExportString and
 * the .json branch of Import / Export.
 *
 * Two element mappings, as in Mathematica:
 *   "RawJSON"  objects <-> Association, arrays <-> List
 *   "JSON"     objects <-> List of Rule (the classic rule-list form)
 * Both map strings to String, true/false/null to True/False/Null, integers to
 * Integer (arbitrary size), numbers with a fraction to Real, and numbers with
 * an exponent but no fraction to the exact Integer / Rational they denote.
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

/* Parse `len` bytes of JSON.  `raw` selects the RawJSON mapping.  Returns an
 * owned expression, or NULL on a syntax error after emitting the matching
 * `<msghead>::json...` messages (msghead is "Import" for ImportString). */
Expr* json_parse(const char* text, size_t len, bool raw, const char* msghead);

/* Serialise `e`.  `raw` selects the RawJSON mapping (rule lists are not
 * objects); `compact` drops all whitespace, otherwise the layout matches
 * Mathematica's (tab indentation, "key":value).  Returns a heap string, or
 * NULL after emitting an Export::json... message when `e` has no JSON form. */
char* json_serialize(const Expr* e, bool raw, bool compact);

/* Registers ImportString / ExportString and wraps Import / Export with the
 * JSON formats.  Must run after imageio_init. */
void json_init(void);

#endif /* MATHILDA_JSON_H */
