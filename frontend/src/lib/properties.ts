// properties.ts — state for the properties sidebar.
//
// Two stores, and they are here rather than in theme.ts or status.ts because
// neither is about a theme or a measurement: one is whether the panel is open,
// the other is a display preference the panel owns.
//
// Plain writables with no persistence, matching darkMode in theme.ts. Nothing in
// this app persists UI state yet, and adding localStorage for one preference
// would make this the only setting that survives a restart -- a surprise, not a
// feature. When persistence arrives it should arrive for all of them at once.

import { writable } from 'svelte/store';

/** Whether the properties sidebar is showing. Toolbar's Sidebar group toggles it. */
export const propertiesOpen = writable(false);

/** Whether a code cell shows its `In[n]` label once it has been evaluated.
 *
 *  Lives here rather than in the toolbar because it is a persistent preference
 *  about how the notebook reads, not an action -- the toolbar is for verbs. */
export const showExecLabels = writable(true);

/** Global UI scale, applied by App.svelte as the root font size.
 *
 *  Lives here because the properties panel offers it, but it is NOT new: the
 *  Cmd+= / Cmd+- / Cmd+0 bindings have driven this value since the toolbar
 *  landed, as a local in App.svelte. Lifting it to a store is what lets the panel
 *  show the same number the keyboard changes -- a second, independent scale would
 *  drift from the first the moment either was used.
 *
 *  The steps the panel offers are exactly representable in binary (0.75, 1, 1.25,
 *  1.5), so a comparison against one is an exact equality rather than an epsilon.
 *  The keyboard still moves in 0.1 increments, which is why the panel highlights a
 *  step only when the value is exactly on it. */
export const uiScale = writable(1.0);

/** The steps the panel offers. Not the keyboard's 0.1: a discrete set is what a
 *  set of buttons can honestly represent. */
export const UI_SCALE_STEPS = [0.75, 1, 1.25, 1.5];
