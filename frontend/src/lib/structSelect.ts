/* structSelect.ts — bottom-up structural selection over bracket/token nesting.
 *
 * Pure, DOM-free, CodeMirror-free: given the token list from mathildaLex.tokenize,
 * build the set of selectable spans, then expand/shrink through them. This is the
 * confirmed frontend-only design — structure comes from brackets and tokens, NOT
 * operator precedence. So in `a + b*c` the ladder from `b` is `b` -> `a + b*c`
 * (there is deliberately no `b*c` step; that would need the kernel parser).
 *
 * Spans are {from, to} in UTF-16 string offsets, matching CodeMirror positions
 * directly. */

import type { Token } from './mathildaLex';

export interface Span { from: number; to: number; }

const OPENERS: Record<string, string> = { '(': ')', '[': ']', '{': '}', '<|': '|>' };
const CLOSERS: Record<string, string> = { ')': '(', ']': '[', '}': '{', '|>': '<|' };

/* "Word-like" tokens are preferred over punctuation when a click lands exactly on
 * a token boundary — double-clicking the edge of `Sin` should select `Sin`, not
 * the `[` beside it. */
const WORD_KINDS = new Set<string>([
  'symbolBuiltin', 'symbolVar', 'number', 'string', 'pattern', 'slot', 'out', 'comment',
]);

function size(s: Span): number { return s.to - s.from; }

function dedup(spans: Span[]): Span[] {
  const seen = new Set<string>();
  const out: Span[] = [];
  for (const s of spans) {
    if (s.to <= s.from) continue;
    const k = s.from + ':' + s.to;
    if (!seen.has(k)) { seen.add(k); out.push(s); }
  }
  return out;
}

/* Produce every selectable span for a document:
 *   - each token (the atoms — the bottom of every ladder);
 *   - for each matched bracket pair, the contents (without delimiters) and the
 *     whole group (with delimiters, plus the head symbol for a call `f[...]`);
 *   - each top-level ;-separated statement, and the whole expression.
 * Brackets are matched as single characters (plus <| |>), which sidesteps the
 * `]]` ambiguity a pure lexer cannot resolve; the only cost is at most one extra
 * nesting level on genuine [[ ]] Part expressions, harmless for selection. */
export function buildSpans(tokens: Token[], src: string): Span[] {
  const spans: Span[] = [];
  const push = (from: number, to: number) => { if (to > from) spans.push({ from, to }); };

  for (const tk of tokens) push(tk.start, tk.end);

  const stack: number[] = [];
  for (let idx = 0; idx < tokens.length; idx++) {
    const tk = tokens[idx];
    if (tk.kind !== 'bracket') continue;
    const text = src.slice(tk.start, tk.end);
    if (OPENERS[text]) { stack.push(idx); continue; }
    if (!CLOSERS[text] || stack.length === 0) continue;
    const openIdx = stack.pop()!;
    const openTk = tokens[openIdx];
    const openText = src.slice(openTk.start, openTk.end);
    push(openTk.end, tk.start);              // contents, delimiters excluded
    let groupFrom = openTk.start;            // whole group, delimiters included
    if (openText === '[' && openIdx > 0) {   // a call/part: include the head symbol
      const prev = tokens[openIdx - 1];
      if (prev.kind === 'symbolBuiltin' || prev.kind === 'symbolVar') groupFrom = prev.start;
    }
    push(groupFrom, tk.end);
  }

  // top-level statements (split on a lone ; at bracket depth 0) + the whole expr
  let depth = 0, segStart = -1, lastEnd = -1;
  for (const tk of tokens) {
    const text = tk.kind === 'bracket' ? src.slice(tk.start, tk.end) : '';
    if (OPENERS[text]) depth++;
    else if (CLOSERS[text] && depth > 0) depth--;
    if (tk.kind === 'semicolon' && depth === 0) {
      if (segStart >= 0) push(segStart, lastEnd);
      segStart = -1;
      continue;
    }
    if (segStart < 0) segStart = tk.start;
    lastEnd = tk.end;
  }
  if (segStart >= 0) push(segStart, lastEnd);
  if (tokens.length) push(tokens[0].start, tokens[tokens.length - 1].end);

  return dedup(spans);
}

/* Pick the atom span under a click/caret `point`. A point strictly inside a token
 * takes it; on a boundary, prefer a word token to the right, then a word to the
 * left, then any token to the right, then to the left. Returns null if there are
 * no tokens at the point (e.g. empty document). */
export function pickAtom(tokens: Token[], point: number): Span | null {
  let interior: Token | null = null;
  let rightWord: Token | null = null, leftWord: Token | null = null;
  let right: Token | null = null, left: Token | null = null;
  for (const tk of tokens) {
    if (tk.start < point && point < tk.end) {
      if (!interior || (tk.end - tk.start) < (interior.end - interior.start)) interior = tk;
    } else if (tk.start === point) {
      if (WORD_KINDS.has(tk.kind)) { if (!rightWord) rightWord = tk; } else if (!right) right = tk;
    } else if (tk.end === point) {
      if (WORD_KINDS.has(tk.kind)) { if (!leftWord) leftWord = tk; } else if (!left) left = tk;
    }
  }
  const tk = interior || rightWord || leftWord || right || left;
  return tk ? { from: tk.start, to: tk.end } : null;
}

/* The smallest span strictly larger than `cur` that still contains it — one step
 * up the ladder. Returns null at the top. */
export function chooseExpand(spans: Span[], cur: Span): Span | null {
  let best: Span | null = null;
  for (const s of spans) {
    if (s.from <= cur.from && s.to >= cur.to && size(s) > size(cur)) {
      if (!best || size(s) < size(best)) best = s;
    }
  }
  return best;
}

/* The full bottom-up ladder from an atom span: [atom, expand, expand, ...]. */
export function ladderFrom(spans: Span[], atom: Span): Span[] {
  const chain: Span[] = [atom];
  let cur = atom;
  for (;;) {
    const next = chooseExpand(spans, cur);
    if (!next) break;
    chain.push(next);
    cur = next;
  }
  return chain;
}

/* Convenience: the ladder for a point, picking the atom first. */
export function ladderAt(tokens: Token[], spans: Span[], point: number): Span[] {
  const atom = pickAtom(tokens, point);
  return atom ? ladderFrom(spans, atom) : [];
}

/* ---- Click gesture state machine (pure, so it is unit-tested) ----
 * A "gesture" is one run of progressive selection anchored on the token first
 * clicked. The FIRST click selects that token; each further click on it climbs
 * one enclosing level, up to the whole expression, then wraps back to the token.
 * A click on a DIFFERENT token starts a new gesture there.
 *
 * It stores only {anchor, level} and recomputes the ladder from the CURRENT
 * spans on every click — so when the span source upgrades mid-gesture (the local
 * bracket/token fallback used for the first click, replaced by the kernel's
 * precedence-aware spans once they arrive) the climb simply uses the better
 * ladder. Anchoring on the atom (not "is the click inside the selection") avoids
 * getting stuck once the whole expression is selected, since then every click is
 * inside it. */
export interface GestureState {
  anchor: Span | null; // the atom this gesture is anchored on
  level: number;       // -1 = none; 0 = the atom; climbs toward the whole expr
}

export function emptyGesture(): GestureState {
  return { anchor: null, level: -1 };
}

export function sameSpan(a: Span | null, b: Span | null): boolean {
  return !!a && !!b && a.from === b.from && a.to === b.to;
}

/* `selection === null` means "nothing to select" (empty space / no tokens). */
export interface GestureResult { state: GestureState; selection: Span | null; }

/* One left click. `atom` is pickAtom() at the click point (may be null); `spans`
 * is the best currently-available span set (kernel if ready, else local). */
export function clickGesture(prev: GestureState, spans: Span[], atom: Span | null): GestureResult {
  if (!atom) return { state: emptyGesture(), selection: null };
  const ladder = ladderFrom(spans, atom);
  if (!ladder.length) return { state: emptyGesture(), selection: null };
  let level: number;
  if (sameSpan(atom, prev.anchor)) {
    // Same token again: climb, wrapping back to the atom after the whole expr.
    level = prev.level < 0 ? 0 : (prev.level >= ladder.length - 1 ? 0 : prev.level + 1);
  } else {
    // New token: select it (the atom) as the start of a fresh gesture.
    level = 0;
  }
  return { state: { anchor: atom, level }, selection: ladder[level] };
}

/* One Alt-Up (expand) / Alt-Down (shrink) step, given the editor's current
 * selection. Re-seeds the gesture from `caret`/selection when it is out of sync
 * with where the user actually is. */
export function stepGesture(
  prev: GestureState, spans: Span[], atom: Span | null,
  cur: Span, grow: boolean,
): GestureResult {
  if (!atom) return { state: prev, selection: null };
  const ladder = ladderFrom(spans, atom);
  if (!ladder.length) return { state: prev, selection: null };
  // Find the current level: the gesture's, if it matches the selection, else the
  // ladder entry equal to (or smallest containing) the current selection.
  let level = prev.level;
  if (level < 0 || level >= ladder.length || !sameSpan(ladder[level], cur)) {
    level = ladder.findIndex((sp) => sameSpan(sp, cur));
    if (level < 0) level = ladder.findIndex((sp) => sp.from <= cur.from && sp.to >= cur.to);
    if (level < 0) level = 0;
  }
  const next = grow ? Math.min(level + 1, ladder.length - 1) : Math.max(level - 1, 0);
  return { state: { anchor: atom, level: next }, selection: ladder[next] };
}
