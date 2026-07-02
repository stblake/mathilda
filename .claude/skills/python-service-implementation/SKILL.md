---
name: python-service-implementation
description: Use when the user asks to implement a new feature, endpoint, service, or background job in a Python backend (FastAPI, Flask, Django, or a standalone worker). Covers the full loop: read the existing code, implement, test, and prepare for review. Trigger on "add an endpoint", "implement a service", "write a worker for X".
---

# python-service-implementation

End-to-end execution of a Python backend change. Opinionated toward FastAPI
+ async + pydantic v2, but adapts to whatever's already in the project.

## When to invoke

- User asks to add / extend / refactor a Python backend feature.
- User references a FastAPI router, a Django view, a Celery task, or similar.
- User says "wire up the database layer", "add the webhook handler", etc.

## Method

1. **Read the project.** Check `pyproject.toml` / `requirements.txt` for the
   framework and ORM in use. Skim `README` and any `CLAUDE.md` / `AGENTS.md`.
   Match their conventions — don't introduce new ones.
2. **Find the natural home.** Which module does this feature belong in?
   Extend existing files before creating new ones. Default layout when
   creating: `api/` (routers), `services/` (logic), `models/` (pydantic +
   ORM), `db/` (session), `core/` (config, logging).
3. **Implement boundary-first.** Define the Pydantic request/response model,
   then the router signature, then the service function, then the data
   access. Each layer gets the types from the layer below it.
4. **Test the behavior.** `pytest` + `pytest-asyncio` + `httpx.AsyncClient`.
   Don't mock the DB when integration-testing a migration path — use a real
   test database.
5. **Verify locally.** Run the test suite and the typechecker; run the
   service briefly against a sample request if the project supports it.
6. **Open the PR** via [`review-pr`](../../commands/review-pr.md) or the
   project's PR command.

## Conventions

- **Async all the way down** when the framework is async. Don't mix
  `requests` into an async handler; use `httpx.AsyncClient`.
- **Pydantic v2.** `BaseModel` + `model_config = ConfigDict(...)`.
- **Types everywhere at boundaries.** Internal helpers can stay inferred
  when obvious.
- **Errors raise at the service layer**; routers map them to HTTPException
  via exception handlers — not inline try/except.
- **Config from `pydantic-settings`**, not hardcoded.

## Skeleton — a FastAPI endpoint

```python
# api/users.py
from fastapi import APIRouter, Depends
from app.schemas.users import UserCreate, UserRead
from app.services.users import UserService

router = APIRouter(prefix="/users", tags=["users"])

@router.post("", response_model=UserRead, status_code=201)
async def create_user(
    payload: UserCreate,
    users: UserService = Depends(),
) -> UserRead:
    return await users.create(payload)
```

```python
# services/users.py
from app.db.session import Session
from app.models.users import User
from app.schemas.users import UserCreate, UserRead

class UserService:
    def __init__(self, db: Session = Depends(get_session)) -> None:
        self.db = db

    async def create(self, payload: UserCreate) -> UserRead:
        user = User(**payload.model_dump())
        self.db.add(user)
        await self.db.commit()
        await self.db.refresh(user)
        return UserRead.model_validate(user)
```

## Guardrails

- Follow [`../../shared/guidelines/security-baseline.md`](../../shared/guidelines/security-baseline.md)
  for secrets, validation, and error handling.
- Follow the engineering principles in root [`../../AGENTS.md`](../../AGENTS.md)
  for minimum-change and boundary validation.
- Don't introduce a new state / query library when the project already uses one.

## What not to do

- No `print()` for debugging. Use `logging`.
- No bare `except Exception:` at router level.
- No silent swallowing — let exceptions reach the handler chain.
- No file-level docstrings that just restate the filename.
