# File I/O

Builtins implemented in `src/readwrite.c`.

## Get
Reads a sequence of Mathilda expressions from a file, evaluates each in order, and returns the value of the last one.
- `Get["filename"]`

**Features**:
- `Protected`.
- Returns `$Failed` if the file cannot be opened.
- Used by the REPL bootstrap to load `src/internal/init.m` (and the rules it pulls in).
- Files conventionally end with `.m`.

## Put
Writes one or more expressions to a file, replacing any prior contents.
- `expr >> "filename"` — shorthand for `Put[expr, "filename"]`.
- `expr >> filename` — bare-word filename, equivalent to `expr >> "filename"`.
- `Put[expr, "filename"]`
- `Put[expr_1, expr_2, ..., expr_n, "filename"]` — writes the expressions one per line.
- `Put["filename"]` — creates an empty file (or truncates an existing one).

**Features**:
- `Protected`.
- The last argument must be a string; it is interpreted as a filename.
- Each `expr_i` is rendered with the standard printer (the same form used at the REPL) and followed by a single `\n`.
- Truncates the file before writing — preserves nothing from a prior `Put`/`PutAppend`.
- Returns `Null` on success and `$Failed` on I/O error.

**Parser notes**:
- `>>` has a low precedence (`30`) — below `Set`/`SetDelayed` (`40`) and above `CompoundExpression` (`10`). `a + b >> "f"` therefore parses as `Put[a + b, "f"]`.
- The bare-word filename accepts identifier characters plus `.`, `/`, `\`, `-`, `_`, `~`, and `$`.

## PutAppend
Like `Put`, but appends to the file rather than truncating it.
- `expr >>> "filename"` — shorthand for `PutAppend[expr, "filename"]`.
- `expr >>> filename` — bare-word filename, equivalent to `expr >>> "filename"`.
- `PutAppend[expr, "filename"]`
- `PutAppend[expr_1, ..., expr_n, "filename"]`

**Features**:
- `Protected`.
- Creates the file if it does not exist; otherwise preserves prior contents and appends new lines.
- Returns `Null` on success and `$Failed` on I/O error.

**Example**:
```
FactorInteger[40320]      >>  "factorizations"
FactorInteger[479001600]  >>> "factorizations"
Get["factorizations"]     (* {{2, 10}, {3, 5}, {5, 2}, {7, 1}, {11, 1}} *)
```
