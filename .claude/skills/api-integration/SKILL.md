---
name: api-integration
description: Use when the user asks to integrate a third-party REST or GraphQL API — fetching data, authenticating, handling rate limits, retrying on failure, and wiring responses into the application.
allowed-tools: Read, Grep, Glob, Bash, Edit, Write
---

# api-integration

End-to-end implementation of a third-party API client: authentication,
request/response modelling, error handling, retries, and rate-limit
compliance. Ships a working `client.py` with tests.

## When to invoke

- User says "integrate the X API", "add a client for Y", "call the Z API".
- User shares an API key or bearer token and asks to start fetching data.
- Existing API client has no retry logic, swallows errors, or lacks tests.
- User needs to paginate through a large API response.

## Method

1. **Read the project.** Does an HTTP client already exist? Match the stack
   (`httpx` if async, `requests` if sync; `gql` for GraphQL).
2. **Model the API surface.** Read the docs URL the user provides (or ask).
   Identify: base URL, auth scheme, key endpoints, pagination strategy, rate
   limits, and error shapes.
3. **Scaffold the client** using the skeleton below. One class per API; one
   method per logical operation.
4. **Implement auth.** Inject API key / OAuth token via a header or param.
   Read credentials from env — never hardcode.
5. **Add retry + back-off.** Use `tenacity` for idiomatic Python retries.
   Respect `Retry-After` headers when present.
6. **Handle errors.** Map HTTP status codes and API error bodies to typed
   exceptions. Never return `None` silently.
7. **Paginate.** Cursor-based pagination in a generator; offset-based with
   a `page` loop. Yield items, not pages.
8. **Write tests** using `respx` (httpx) or `responses` (requests) to mock
   the network. Cover: happy path, 429 retry, auth failure, pagination end.
9. **Type everything.** Pydantic models for response bodies.

## Skeleton

```python
# client.py
from __future__ import annotations

import os
from typing import Any, Generator

import httpx
from pydantic import BaseModel
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

BASE_URL = "https://api.example.com/v1"


class APIError(Exception):
    def __init__(self, status: int, detail: str) -> None:
        super().__init__(f"HTTP {status}: {detail}")
        self.status = status


class ExampleItem(BaseModel):
    id: str
    name: str
    # ... add fields from API schema


class ExampleAPIClient:
    def __init__(self, api_key: str | None = None) -> None:
        self._api_key = api_key or os.environ["EXAMPLE_API_KEY"]
        self._client = httpx.AsyncClient(
            base_url=BASE_URL,
            headers={"Authorization": f"Bearer {self._api_key}"},
            timeout=30.0,
        )

    async def __aenter__(self) -> "ExampleAPIClient":
        return self

    async def __aexit__(self, *args: Any) -> None:
        await self._client.aclose()

    @retry(
        retry=retry_if_exception_type(httpx.HTTPStatusError),
        wait=wait_exponential(multiplier=1, min=1, max=60),
        stop=stop_after_attempt(4),
    )
    async def _get(self, path: str, **params: Any) -> dict[str, Any]:
        r = await self._client.get(path, params=params)
        if r.status_code == 429:
            raise httpx.HTTPStatusError("rate limited", request=r.request, response=r)
        if not r.is_success:
            raise APIError(r.status_code, r.text)
        return r.json()

    async def list_items(self, *, page_size: int = 100) -> Generator[ExampleItem, None, None]:
        cursor: str | None = None
        while True:
            params: dict[str, Any] = {"limit": page_size}
            if cursor:
                params["cursor"] = cursor
            data = await self._get("/items", **params)
            for raw in data["items"]:
                yield ExampleItem(**raw)
            cursor = data.get("next_cursor")
            if not cursor:
                break
```

```python
# test_client.py
import pytest
import respx
import httpx
from .client import ExampleAPIClient, ExampleItem

@pytest.mark.asyncio
@respx.mock
async def test_list_items_happy_path():
    respx.get("https://api.example.com/v1/items").mock(
        return_value=httpx.Response(200, json={
            "items": [{"id": "1", "name": "foo"}],
            "next_cursor": None,
        })
    )
    async with ExampleAPIClient(api_key="test-key") as client:
        items = [item async for item in client.list_items()]
    assert len(items) == 1
    assert isinstance(items[0], ExampleItem)

@pytest.mark.asyncio
@respx.mock
async def test_retries_on_429():
    route = respx.get("https://api.example.com/v1/items")
    route.side_effect = [
        httpx.Response(429, json={"error": "rate_limited"}),
        httpx.Response(200, json={"items": [], "next_cursor": None}),
    ]
    async with ExampleAPIClient(api_key="test-key") as client:
        items = [item async for item in client.list_items()]
    assert items == []
```

## Guardrails

- Follow [`../../shared/guidelines/security-baseline.md`](../../shared/guidelines/security-baseline.md):
  credentials from env, never hardcoded.
- Never log raw API responses that may contain PII or secrets.
- Set explicit timeouts — never leave `timeout=None`.
- Pin the API version in the base URL when the provider offers versioning.

## Anti-patterns

- `requests.get(url)` without error handling — always check status.
- Catching `Exception` broadly and returning `None` — raises silently fail.
- Building pagination with `while True` and no exit condition guard.
- Storing tokens in source files or config committed to git.

## Extending

- GraphQL: swap `_get` for a `_query(document, variables)` method using `gql`.
- OAuth 2.0: add a token-refresh flow with a `_refresh_token` method.
- Webhook ingestion: add a `verify_signature(payload, header)` helper.
- Caching: wrap reads with `functools.lru_cache` or a Redis layer for quota-heavy endpoints.
