/* mathildaLang.ts — CodeMirror 6 syntax highlighting for Mathilda source.
 *
 * Wraps the pure lexer in mathildaLex.ts as a StreamLanguage (per-line tokeniser
 * that carries state across newlines — exactly what nested comments and
 * multi-line strings need) and pairs it with a HighlightStyle whose colours are
 * CSS variables, so the single style themes itself for light/dark the same way
 * the rest of the UI does (the --cm-* vars are defined in App.svelte).
 *
 * StreamLanguage is used in preference to a Lezer grammar: Mathilda has no Lezer
 * grammar, and highlighting needs only lexemes, not a parse tree. */

import { StreamLanguage, LanguageSupport, HighlightStyle, syntaxHighlighting } from '@codemirror/language';
import { tags as t, type Tag } from '@lezer/highlight';
import { scanToken, startLexState } from './mathildaLex';
import type { LexState } from './mathildaLex';

/* Token-kind name (returned by the stream tokeniser) → highlight tag. The names
 * match TokenKind in mathildaLex.ts; 'ws' returns null (no tag) instead. */
const tokenTable: Record<string, Tag> = {
  comment: t.blockComment,
  string: t.string,
  number: t.number,
  symbolBuiltin: t.function(t.variableName),
  symbolVar: t.variableName,
  pattern: t.special(t.variableName),
  slot: t.meta,
  out: t.atom,
  operator: t.operator,
  semicolon: t.operator,
  bracket: t.bracket,
  comma: t.separator,
  error: t.invalid,
};

const mathildaStream = StreamLanguage.define<LexState>({
  name: 'mathilda',
  startState: () => startLexState(),
  copyState: (s) => ({ comment: s.comment, inString: s.inString }),
  token(stream, state) {
    const line = stream.string;
    const start = stream.pos;
    const r = scanToken(line, start, line.length, state);
    stream.pos = r.end > start ? r.end : start + 1; // never stall
    state.comment = r.state.comment;
    state.inString = r.state.inString;
    return r.kind === 'ws' ? null : r.kind;
  },
  tokenTable,
});

/** The language extension to add to an EditorView. */
export const mathildaLanguage = new LanguageSupport(mathildaStream);

/** One HighlightStyle; colours resolve from --cm-* CSS vars (App.svelte). */
export const mathildaHighlight = HighlightStyle.define([
  { tag: t.blockComment, color: 'var(--cm-comment)', fontStyle: 'italic' },
  { tag: t.string, color: 'var(--cm-string)' },
  { tag: t.number, color: 'var(--cm-number)' },
  { tag: t.function(t.variableName), color: 'var(--cm-builtin)' },
  { tag: t.variableName, color: 'var(--cm-symbol)' },
  { tag: t.special(t.variableName), color: 'var(--cm-pattern)' },
  { tag: t.meta, color: 'var(--cm-slot)' },
  { tag: t.atom, color: 'var(--cm-out)' },
  { tag: t.operator, color: 'var(--cm-operator)' },
  { tag: t.bracket, color: 'var(--cm-bracket)' },
  { tag: t.separator, color: 'var(--cm-bracket)' },
  { tag: t.invalid, color: 'var(--cm-error)' },
]);

/** Drop-in extension bundle: language + its highlight style. */
export const mathildaHighlightExtension = [mathildaLanguage, syntaxHighlighting(mathildaHighlight)];
