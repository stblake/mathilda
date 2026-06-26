---
name: react-feature-delivery
description: Use when the user asks to add a component, screen, feature, or user interaction in a React + TypeScript frontend (also applicable to Vue when the project uses Vue). Covers implementation, state and data wiring, accessibility, and tests. Trigger on "add a button that…", "build the X screen", "wire this form up".
---

# react-feature-delivery

Full-loop delivery of a React + TypeScript UI change. Adapts to the project's
existing styling, state, and query layer — does not introduce new ones.

## When to invoke

- User asks to add, change, or fix UI in a React or Vue project.
- User describes user-facing behavior ("when the user clicks X, Y happens").
- User requests a new route / page / modal / form.

## Method

1. **Read the project.** Look at `package.json` — what's the state story
   (Zustand, Redux Toolkit, Jotai, context), query layer (TanStack Query,
   SWR, Apollo), styling (Tailwind, CSS Modules, styled-components)?
   Match all three.
2. **Find the natural home.** Which feature folder does this live in?
   Extend existing components before creating new ones. When creating:
   `components/` (presentational), `features/<feature>/` (business logic),
   `hooks/` (reusable hooks).
3. **Implement props-first.** Declare the `interface FooProps` before the
   component body. Named export, function component, no `React.FC`.
4. **Wire state and data last.** Local `useState` first. Reach for a store
   only if state is shared across unrelated subtrees. Server state through
   the query layer.
5. **Accessibility on the first pass.** Right element, right label, right
   focus order. Don't retrofit later.
6. **Test behavior.** Vitest / Jest + React Testing Library. Query by role
   or text. Mock at the network boundary (MSW), not at the query layer.
7. **Verify in the browser.** Run the dev server and drive the feature by
   hand — type checks alone do not prove the feature works.

## Skeleton — a typed component

```tsx
// features/users/UserCard.tsx
interface UserCardProps {
  user: { id: string; name: string; email: string };
  onEdit?: (id: string) => void;
}

export function UserCard({ user, onEdit }: UserCardProps) {
  return (
    <article className="rounded-lg border p-4">
      <h3 className="font-semibold">{user.name}</h3>
      <p className="text-sm text-muted-foreground">{user.email}</p>
      {onEdit && (
        <button type="button" onClick={() => onEdit(user.id)} aria-label={`Edit ${user.name}`}>
          Edit
        </button>
      )}
    </article>
  );
}
```

## Skeleton — a data hook

```ts
// features/users/useUser.ts
import { useQuery } from '@tanstack/react-query';

export function useUser(id: string) {
  return useQuery({
    queryKey: ['user', id],
    queryFn: () => fetch(`/api/users/${id}`).then(r => r.json()),
  });
}
```

## Guardrails

- Follow [`../../shared/guidelines/security-baseline.md`](../../shared/guidelines/security-baseline.md)
  — no `dangerouslySetInnerHTML` with untrusted input, no secrets in the bundle.
- Follow the engineering principles in the library root [`../../AGENTS.md`](../../AGENTS.md)
  — extend before create, no speculative abstractions.

## What not to do

- No new `useEffect` for data fetching. Use the query layer.
- No `any` types. `unknown` + narrowing instead.
- No `// eslint-disable` without a reason comment.
- No default exports for new components.
- No committed `console.log`.

## Vue variant

For Vue projects, same method, different primitives:

- `<script setup lang="ts">` + `defineProps<Props>()`.
- Composables in `composables/` replace custom hooks.
- Pinia for shared state; `ref` / `reactive` for local.
