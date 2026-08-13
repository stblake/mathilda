// theme.ts — light/dark choice.
//
// Hoisted out of App.svelte because the toolbar's Notebook group and the
// properties panel both need to read and toggle it, and prop-drilling a store
// through a component that only forwards it is noise.
//
// App.svelte still owns the side effect of putting the class on <html>: that is
// a DOM write and belongs in a component's reactive statement, not here.

import { writable } from 'svelte/store';

/** The canvas is dark by default. */
export const darkMode = writable(true);
