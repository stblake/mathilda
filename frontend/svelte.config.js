import { vitePreprocess } from '@sveltejs/vite-plugin-svelte';

/** @type {import("@sveltejs/vite-plugin-svelte").SvelteConfig} */
export default {
  preprocess: vitePreprocess(),
  compilerOptions: {
    // Use Svelte 4-compatible (legacy) mode so on:click / on:event /
    // createEventDispatcher / $: reactive labels all work as expected.
    runes: false,
  },
};
