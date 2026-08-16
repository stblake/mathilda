// snippets.ts — the Insert group's template palette.
//
// Each entry is a CodeMirror snippet template: `${name}` marks a field, and Tab
// walks them after insertion. `snippet()` from @codemirror/autocomplete does the
// work; this module is the templates plus the two functions the toolbar needs.
//
// The templates are the risk, not the machinery. A template that inserts an
// unbalanced bracket is a broken cell every time it is used, and reading one is a
// poor way to notice -- `{{${a}, ${b}}, {${c}, ${d}}}` has to be counted. So
// `expandedText` exists, `npm run check:snippets` feeds each expansion through the
// real Mathilda binary, and a template whose expansion does not PARSE fails the
// check. That is a stronger guarantee than balanced brackets: it is the language's
// own parser agreeing that what the button inserts is a valid expression.

import { snippet } from '@codemirror/autocomplete';
import type { EditorView } from '@codemirror/view';

export type SnippetDef = {
  /** Menu label. */
  label: string;
  /** CodeMirror snippet template; `${name}` is a tab field. */
  template: string;
};

/* Ordered by how often the thing is reached for, not alphabetically: a palette is
 * read top-down and the first two entries should be the ones that save the most
 * typing. Every expansion is checked against Mathilda's parser. */
export const SNIPPETS: SnippetDef[] = [
  { label: 'Table',       template: 'Table[${expr}, {${i}, 1, ${n}}]' },
  { label: 'Matrix 2×2',  template: '{{${a}, ${b}}, {${c}, ${d}}}' },
  { label: 'Sum',         template: 'Sum[${expr}, {${i}, 1, ${n}}]' },
  { label: 'Product',     template: 'Product[${expr}, {${i}, 1, ${n}}]' },
  { label: 'Integrate',   template: 'Integrate[${expr}, ${x}]' },
  { label: 'NIntegrate',  template: 'NIntegrate[${expr}, {${x}, 0, 1}]' },
  { label: 'Derivative',  template: 'D[${expr}, ${x}]' },
  { label: 'Limit',       template: 'Limit[${expr}, ${x} -> 0]' },
  { label: 'Solve',       template: 'Solve[${lhs} == ${rhs}, ${x}]' },
  { label: 'Plot',        template: 'Plot[${expr}, {${x}, -5, 5}]' },
  { label: 'Module',      template: 'Module[{${vars}}, ${body}]' },
  { label: 'Definition',  template: '${f}[${x}_] := ${body}' },
  { label: 'Replace all', template: '${expr} /. ${lhs} -> ${rhs}' },
  { label: 'Map',         template: 'Map[${f}, ${list}]' },
];

/** The text a template becomes once every field is left at its placeholder name.
 *
 *  What the user sees the instant the snippet lands, so it doubles as the menu's
 *  hint and as what the parser check parses. `${expr}` becomes `expr`. */
export function expandedText(template: string): string {
  return template.replace(/\$\{([^}]*)\}/g, '$1');
}

/** Insert `template` at the caret, replacing any selection.
 *
 *  Declines when there is no editor rather than guessing a target: the Insert
 *  group is only shown for a code cell, but the active cell can be deleted between
 *  the menu opening and an item being chosen. */
export function insertSnippet(view: EditorView | null, template: string): boolean {
  if (!view) return false;
  const { from, to } = view.state.selection.main;
  /* snippet() returns an apply function shaped for the completion system, which
     passes the completion it came from; there is no completion here, and it is
     unused for a plain insertion. */
  snippet(template)(view, null as never, from, to);
  view.focus();
  return true;
}
