#ifndef LIST_H
#define LIST_H

/* Umbrella header for the list builtins, split across src/list/.
 *
 * External callers continue to `#include "list.h"` unchanged; this file
 * re-exports every per-module public header plus the list_init() entry point
 * that core_init() calls. */

#include "expr.h"

#include "table.h"
#include "range.h"
#include "array.h"
#include "constant_array.h"
#include "take_drop.h"
#include "flatten.h"
#include "partition.h"
#include "rotate.h"
#include "reverse.h"
#include "transpose.h"
#include "setops.h"
#include "split.h"
#include "rescale.h"
#include "pad.h"
#include "total.h"
#include "accumulate.h"
#include "differences.h"
#include "ratios.h"
#include "listpredicates.h"
#include "matrixq.h"
#include "minmax.h"
#include "join.h"

void list_init(void);

#endif // LIST_H
