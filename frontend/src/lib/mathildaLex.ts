/* mathildaLex.ts — a pure, DOM-free lexer for Mathilda (Wolfram-style) source.
 *
 * This is NOT a parser: it recognises lexemes only (no precedence, no tree). It
 * is the single source of truth for both features that need to read code:
 *   • syntax highlighting (mathildaLang.ts wraps it in a CodeMirror StreamLanguage)
 *   • bottom-up structural selection (structSelect.ts builds spans from the tokens)
 *
 * The rules mirror the kernel's inline lexer in src/parse.c (parse_symbol,
 * parse_number, parse_string, get_operator, the pattern/slot cases in
 * parse_primary, and the nested (* *) comment handling in skip_whitespace) so the
 * notebook colours and selects exactly what the kernel would read. It deliberately
 * does not import anything: node can run it directly for tests, and it carries no
 * CodeMirror dependency.
 *
 * Offsets are UTF-16 string indices (what both JS strings and CodeMirror use), so
 * no byte/codepoint translation is ever required — this never crosses to the C
 * kernel.
 */

export type TokenKind =
  | 'comment' | 'string' | 'number'
  | 'symbolBuiltin' | 'symbolVar' | 'pattern'
  | 'slot' | 'out' | 'operator' | 'semicolon'
  | 'bracket' | 'comma' | 'error' | 'ws';

export interface Token { start: number; end: number; kind: TokenKind; }

/* Carried across lines so a (* *) comment or a "..." string can span newlines.
 * `comment` is the nesting depth (0 = not in a comment). */
export interface LexState { comment: number; inString: boolean; }
export function startLexState(): LexState { return { comment: 0, inString: false }; }

export interface ScanResult { kind: TokenKind; end: number; state: LexState; }

const isDigit = (c: string) => c >= '0' && c <= '9';
const isSymStart = (c: string) =>
  (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c === '$';
const isSymChar = (c: string) => isSymStart(c) || isDigit(c) || c === '`';
const isSpace = (c: string) =>
  c === ' ' || c === '\t' || c === '\r' || c === '\n' || c === '\f' || c === '\v';

/* Operator characters. We tokenise a MAXIMAL RUN of these as one `operator`
 * token rather than matching the ~80-entry operator table from get_operator.
 * Rationale: the exact operator boundary never affects highlighting (every
 * operator gets one colour) nor selection (operators are not structural
 * delimiters), while a hand-maintained partial table risks mis-splitting a real
 * operator. `,` and `;` are handled separately because selection needs them. */
const OP_CHARS = '+-*/^.@!?:=<>~&|\\';

/* Consume a numeric literal starting at i: integer/real/scientific, plus the
 * Mathematica precision/accuracy backtick suffix (`n, ``n) and scaled form *^. */
function scanNumber(src: string, i: number, end: number): number {
  let j = i;
  while (j < end && isDigit(src[j])) j++;
  // decimal point, but not the Repeated operator `..`
  if (src[j] === '.' && src[j + 1] !== '.') { j++; while (j < end && isDigit(src[j])) j++; }
  // exponent
  if ((src[j] === 'e' || src[j] === 'E') &&
      (isDigit(src[j + 1]) ||
       ((src[j + 1] === '+' || src[j + 1] === '-') && isDigit(src[j + 2])))) {
    j++; if (src[j] === '+' || src[j] === '-') j++;
    while (j < end && isDigit(src[j])) j++;
  }
  // precision (`n) or accuracy (``n) backtick suffix
  if (src[j] === '`') {
    j++; if (src[j] === '`') j++;
    while (j < end && isDigit(src[j])) j++;
    if (src[j] === '.' && isDigit(src[j + 1])) { j++; while (j < end && isDigit(src[j])) j++; }
  }
  // scaled scientific *^
  if (src[j] === '*' && src[j + 1] === '^') {
    j += 2; if (src[j] === '+' || src[j] === '-') j++;
    while (j < end && isDigit(src[j])) j++;
  }
  return j;
}

/* Consume a "..." body starting at the opening quote; returns index past the
 * closing quote, or `end` if unterminated. */
function scanStringBody(src: string, j: number, end: number): number {
  j++; // opening quote
  while (j < end) {
    if (src[j] === '\\' && j + 1 < end) { j += 2; continue; }
    if (src[j] === '"') { j++; break; }
    j++;
  }
  return j;
}

/* Consume a blank/pattern tail starting at the first `_`:
 * _ / __ / ___ , an optional head (_h, x__h), and an optional trailing `.`
 * (Optional, as in x_.). */
function scanBlankTail(src: string, k: number, end: number): number {
  let j = k, u = 0;
  while (j < end && src[j] === '_' && u < 3) { j++; u++; }
  if (j < end && isSymStart(src[j])) { j++; while (j < end && isSymChar(src[j])) j++; }
  if (src[j] === '.' && src[j + 1] !== '.') j++;
  return j;
}

/* Capitalisation heuristic (confirmed design choice): a name whose base (after
 * the last context backtick) begins with an uppercase letter or `$` is treated
 * as a builtin — the Wolfram convention (Sin, Integrate, $Line) — otherwise a
 * user variable. */
function symbolKind(name: string): TokenKind {
  const bt = name.lastIndexOf('`');
  const base = bt >= 0 ? name.slice(bt + 1) : name;
  const f = base.charCodeAt(0);
  const isBuiltin = (f >= 65 && f <= 90) || base[0] === '$';
  return isBuiltin ? 'symbolBuiltin' : 'symbolVar';
}

/* Scan exactly one token from src starting at i, not reading past `end` (which is
 * the end of the current line for the StreamLanguage, or src.length for a whole-
 * document tokenize). Returns the token kind, the end offset, and the updated
 * cross-line state. Always advances at least one position for a non-empty input,
 * so callers cannot loop forever. */
export function scanToken(src: string, i: number, end: number, state: LexState): ScanResult {
  // continuation: inside a nested comment
  if (state.comment > 0) {
    let depth = state.comment;
    while (i < end) {
      if (src[i] === '(' && src[i + 1] === '*') { depth++; i += 2; }
      else if (src[i] === '*' && src[i + 1] === ')') { depth--; i += 2; if (depth === 0) break; }
      else i++;
    }
    return { kind: 'comment', end: i, state: { comment: depth, inString: false } };
  }
  // continuation: inside a string
  if (state.inString) {
    while (i < end) {
      if (src[i] === '\\' && i + 1 < end) { i += 2; continue; }
      if (src[i] === '"') { i++; return { kind: 'string', end: i, state: { comment: 0, inString: false } }; }
      i++;
    }
    return { kind: 'string', end: i, state: { comment: 0, inString: true } };
  }

  const c = src[i];

  // comment open
  if (c === '(' && src[i + 1] === '*') {
    let depth = 1; i += 2;
    while (i < end) {
      if (src[i] === '(' && src[i + 1] === '*') { depth++; i += 2; }
      else if (src[i] === '*' && src[i + 1] === ')') { depth--; i += 2; if (depth === 0) break; }
      else i++;
    }
    return { kind: 'comment', end: i, state: { comment: depth, inString: false } };
  }
  // string open
  if (c === '"') {
    const j = scanStringBody(src, i, end);
    const closed = j <= end && src[j - 1] === '"' && j > i + 1;
    return { kind: 'string', end: j, state: { comment: 0, inString: !closed } };
  }
  // whitespace
  if (isSpace(c)) { let j = i + 1; while (j < end && isSpace(src[j])) j++; return { kind: 'ws', end: j, state }; }
  // number
  if (isDigit(c) || (c === '.' && isDigit(src[i + 1]))) {
    return { kind: 'number', end: scanNumber(src, i, end), state };
  }
  // slot: #, ##, #n, ##n, #name, #"name"
  if (c === '#') {
    let j = i + 1;
    if (src[j] === '#') { j++; while (j < end && isDigit(src[j])) j++; }
    else if (src[j] === '"') { j = scanStringBody(src, j, end); }
    else if (isDigit(src[j])) { while (j < end && isDigit(src[j])) j++; }
    else if (isSymStart(src[j])) { while (j < end && isSymChar(src[j])) j++; }
    return { kind: 'slot', end: j, state };
  }
  // out: %, %%, %%%, %n
  if (c === '%') {
    let j = i + 1;
    if (src[j] === '%') { j++; while (j < end && src[j] === '%') j++; }
    else while (j < end && isDigit(src[j])) j++;
    return { kind: 'out', end: j, state };
  }
  // leading-blank pattern: _, __, ___ (with optional head / trailing .)
  if (c === '_') return { kind: 'pattern', end: scanBlankTail(src, i, end), state };
  // symbol, or a named pattern (x_, x__h, x_.)
  if (isSymStart(c)) {
    let j = i + 1; while (j < end && isSymChar(src[j])) j++;
    if (src[j] === '_') return { kind: 'pattern', end: scanBlankTail(src, j, end), state };
    return { kind: symbolKind(src.slice(i, j)), end: j, state };
  }
  // association delimiters and ordinary brackets (single chars; see structSelect)
  if (c === '<' && src[i + 1] === '|') return { kind: 'bracket', end: i + 2, state };
  if (c === '|' && src[i + 1] === '>') return { kind: 'bracket', end: i + 2, state };
  if (c === '(' || c === ')' || c === '[' || c === ']' || c === '{' || c === '}')
    return { kind: 'bracket', end: i + 1, state };
  // comma
  if (c === ',') return { kind: 'comma', end: i + 1, state };
  // semicolon: ;; is the Span operator; a lone ; separates statements
  if (c === ';') {
    if (src[i + 1] === ';') return { kind: 'operator', end: i + 2, state };
    return { kind: 'semicolon', end: i + 1, state };
  }
  // maximal operator run
  if (OP_CHARS.indexOf(c) >= 0) {
    let j = i + 1; while (j < end && OP_CHARS.indexOf(src[j]) >= 0) j++;
    return { kind: 'operator', end: j, state };
  }
  // anything else: a single-char error token keeps the lexer total and robust
  return { kind: 'error', end: i + 1, state };
}

/* Tokenise a whole document into non-whitespace tokens, tracking cross-line
 * state so multi-line comments and strings become single tokens. Used by
 * structural selection and by the tests. */
export function tokenize(src: string): Token[] {
  const tokens: Token[] = [];
  const state = startLexState();
  const end = src.length;
  let i = 0;
  while (i < end) {
    const r = scanToken(src, i, end, state);
    const next = r.end > i ? r.end : i + 1; // guarantee forward progress
    if (r.kind !== 'ws') tokens.push({ start: i, end: next, kind: r.kind });
    state.comment = r.state.comment;
    state.inString = r.state.inString;
    i = next;
  }
  return tokens;
}
