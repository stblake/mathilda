/*
 * schemes.ts — syntax-highlighting colour schemes for notebook code cells.
 *
 * The highlighter (mathildaLang.ts) paints every token from a `--cm-*` CSS
 * variable, so a scheme is simply a set of those variables. `applyColorScheme`
 * sets them INLINE on <html>, which beats the stylesheet `:root` / `html.light`
 * rules in App.svelte — so a scheme overrides both light and dark. The `default`
 * scheme defines none, clearing the inline values so App.svelte's adaptive
 * Catppuccin light/dark palette shows through again.
 *
 * The nine fixed palettes are theme-independent (chosen to read on either
 * background); `default` is the one that adapts to light/dark. Two of the fixed
 * ones (Solarized Light, GitHub Light) are tuned for a light canvas.
 */

/** The token roles, matching the `--cm-<key>` variables in App.svelte and the
 *  tags wired up in mathildaLang.ts. */
export const CM_KEYS = [
  'comment', 'string', 'number', 'builtin', 'symbol',
  'pattern', 'slot', 'out', 'operator', 'bracket', 'error',
] as const;
export type CmKey = (typeof CM_KEYS)[number];

export interface ColorScheme {
  id: string;
  label: string;
  /** Token colours. Omitted for `default`, which keeps the adaptive palette.
   *  Values may be hex (fixed palettes), a `var(--…)` reference (adaptive
   *  schemes like grayscale, which resolve against light/dark), or `inherit`
   *  (the `off` scheme, which renders every token in the plain text colour). */
  cm?: Partial<Record<CmKey, string>>;
  /** Render builtin names in bold (used by grayscale, where weight — not hue —
   *  distinguishes them). Drives --cm-builtin-weight; see mathildaLang.ts. */
  builtinBold?: boolean;
}

export const COLOR_SCHEMES: ColorScheme[] = [
  { id: 'default', label: 'Default (adaptive)' },
  {
    id: 'dracula', label: 'Dracula',
    cm: { comment: '#6272a4', string: '#f1fa8c', number: '#bd93f9', builtin: '#8be9fd',
          symbol: '#f8f8f2', pattern: '#ff79c6', slot: '#ffb86c', out: '#ffb86c',
          operator: '#ff79c6', bracket: '#f8f8f2', error: '#ff5555' },
  },
  {
    id: 'nord', label: 'Nord',
    cm: { comment: '#616e88', string: '#a3be8c', number: '#b48ead', builtin: '#88c0d0',
          symbol: '#d8dee9', pattern: '#81a1c1', slot: '#ebcb8b', out: '#ebcb8b',
          operator: '#81a1c1', bracket: '#8891a5', error: '#bf616a' },
  },
  {
    id: 'monokai', label: 'Monokai',
    cm: { comment: '#88846f', string: '#e6db74', number: '#ae81ff', builtin: '#66d9ef',
          symbol: '#f8f8f2', pattern: '#f92672', slot: '#fd971f', out: '#fd971f',
          operator: '#f92672', bracket: '#f8f8f2', error: '#f92672' },
  },
  {
    id: 'solarized-dark', label: 'Solarized Dark',
    cm: { comment: '#586e75', string: '#2aa198', number: '#d33682', builtin: '#268bd2',
          symbol: '#93a1a1', pattern: '#cb4b16', slot: '#b58900', out: '#b58900',
          operator: '#859900', bracket: '#839496', error: '#dc322f' },
  },
  {
    id: 'solarized-light', label: 'Solarized Light',
    cm: { comment: '#93a1a1', string: '#2aa198', number: '#d33682', builtin: '#268bd2',
          symbol: '#586e75', pattern: '#cb4b16', slot: '#b58900', out: '#b58900',
          operator: '#859900', bracket: '#657b83', error: '#dc322f' },
  },
  {
    id: 'gruvbox', label: 'Gruvbox',
    cm: { comment: '#928374', string: '#b8bb26', number: '#d3869b', builtin: '#83a598',
          symbol: '#ebdbb2', pattern: '#fb4934', slot: '#fabd2f', out: '#fabd2f',
          operator: '#fe8019', bracket: '#a89984', error: '#fb4934' },
  },
  {
    id: 'one-dark', label: 'One Dark',
    cm: { comment: '#5c6370', string: '#98c379', number: '#d19a66', builtin: '#61afef',
          symbol: '#abb2bf', pattern: '#c678dd', slot: '#e5c07b', out: '#e5c07b',
          operator: '#56b6c2', bracket: '#abb2bf', error: '#e06c75' },
  },
  {
    id: 'tokyo-night', label: 'Tokyo Night',
    cm: { comment: '#565f89', string: '#9ece6a', number: '#ff9e64', builtin: '#7aa2f7',
          symbol: '#c0caf5', pattern: '#bb9af7', slot: '#e0af68', out: '#e0af68',
          operator: '#89ddff', bracket: '#a9b1d6', error: '#f7768e' },
  },
  {
    id: 'github-light', label: 'GitHub Light',
    cm: { comment: '#6e7781', string: '#0a3069', number: '#0550ae', builtin: '#8250df',
          symbol: '#24292f', pattern: '#cf222e', slot: '#953800', out: '#953800',
          operator: '#cf222e', bracket: '#24292f', error: '#cf222e' },
  },
  {
    // A 1970s green-phosphor CRT look: everything in shades of green, builtins
    // brightest AND bold. Adaptive like grayscale -- the colours are var()
    // references into the --cm-grn-* ramp (App.svelte), which carries bright
    // greens for dark and deep greens for light. Error stays red: a
    // monochrome-green scheme has no green left to signal "wrong" with.
    id: 'seventies', label: '70s (Green Phosphor)',
    builtinBold: true,
    cm: { comment: 'var(--cm-grn-comment)', string: 'var(--cm-grn-string)',
          number: 'var(--cm-grn-number)', builtin: 'var(--cm-grn-builtin)',
          symbol: 'var(--cm-grn-symbol)', pattern: 'var(--cm-grn-pattern)',
          slot: 'var(--cm-grn-slot)', out: 'var(--cm-grn-slot)',
          operator: 'var(--cm-grn-operator)', bracket: 'var(--cm-grn-bracket)',
          error: 'var(--cm-grn-error)' },
  },
  {
    // An 80s amber-terminal look: browns through golds to bright orange.
    id: 'eighties', label: '80s (Brown/Orange)',
    cm: { comment: '#8a5a2b', string: '#d98f00', number: '#ff7a18', builtin: '#ffa033',
          symbol: '#d8a05a', pattern: '#e8622e', slot: '#f0a000', out: '#f0a000',
          operator: '#c9601a', bracket: '#a86a3a', error: '#ff3b1f' },
  },
  {
    // Adaptive greyscale: hue carries no meaning, so builtins are set apart by
    // WEIGHT (bold). The values are var() references into the --cm-gray-* ramp,
    // which App.svelte defines with light and dark variants — hence adaptive.
    id: 'grayscale', label: 'Grayscale',
    builtinBold: true,
    cm: { comment: 'var(--cm-gray-faint)', string: 'var(--cm-gray-mid)',
          number: 'var(--cm-gray-mid)', builtin: 'var(--cm-gray-fg)',
          symbol: 'var(--cm-gray-fg)', pattern: 'var(--cm-gray-dim)',
          slot: 'var(--cm-gray-dim)', out: 'var(--cm-gray-dim)',
          operator: 'var(--cm-gray-strong)', bracket: 'var(--cm-gray-strong)',
          error: 'var(--cm-gray-strong)' },
  },
  {
    // Highlighting off: every token renders in the editor's plain text colour
    // (`inherit`), which follows light/dark on its own.
    id: 'off', label: 'Off (none)',
    cm: { comment: 'inherit', string: 'inherit', number: 'inherit', builtin: 'inherit',
          symbol: 'inherit', pattern: 'inherit', slot: 'inherit', out: 'inherit',
          operator: 'inherit', bracket: 'inherit', error: 'inherit' },
  },
];

/** Apply a scheme by setting (or clearing, for `default`) the `--cm-*` variables
 *  inline on <html>. Inline values override the stylesheet rules, so this wins in
 *  both light and dark; `default` removes them to restore the adaptive palette. */
export function applyColorScheme(id: string): void {
  if (typeof document === 'undefined') return;
  const root = document.documentElement;
  const scheme = COLOR_SCHEMES.find((s) => s.id === id);
  for (const key of CM_KEYS) {
    const color = scheme?.cm?.[key];
    if (color) root.style.setProperty(`--cm-${key}`, color);
    else root.style.removeProperty(`--cm-${key}`);
  }
  /* Builtin weight: bold only for schemes that ask for it (grayscale). Cleared
     otherwise so the highlighter's `var(--cm-builtin-weight, normal)` falls back. */
  if (scheme?.builtinBold) root.style.setProperty('--cm-builtin-weight', 'bold');
  else root.style.removeProperty('--cm-builtin-weight');
}
