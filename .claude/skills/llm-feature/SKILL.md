---
name: llm-feature
description: Use when the user asks to build, integrate, or improve an AI/LLM feature — calling Claude, GPT, or another model API; adding prompt templates; structured output extraction; tool use; streaming; or evaluating model outputs.
allowed-tools: Read, Grep, Glob, Bash, Edit, Write
---

# llm-feature

End-to-end implementation of an LLM-powered feature: prompt engineering,
API wiring, structured output, streaming, tool use, and evaluation. Ships
production patterns — not toy examples.

## When to invoke

- User asks to "add an AI feature", "call Claude/GPT", "use an LLM to…".
- User wants to extract structured data from free-form text using a model.
- User needs to add tool use, function calling, or multi-step agent logic.
- User wants to stream model output to the UI.
- User wants to evaluate or test a prompt systematically.

## Method

1. **Understand the task.** What should the LLM do? What are the inputs?
   What's the expected output shape? What happens when it fails?
2. **Pick the right model.** Fast/cheap (haiku, gpt-4o-mini) for high-volume
   classification; capable (sonnet/opus, gpt-4o) for complex reasoning.
3. **Write the prompt.** See skeleton below. System prompt is the persona
   and constraints; user prompt is the specific task. Separate them.
4. **Add structured output.** Use Pydantic + `instructor` (Python) or
   `response_format` (OpenAI) to get typed, validated responses — never
   parse JSON from raw text manually.
5. **Handle failures.** Models hallucinate, time out, and hit rate limits.
   Validate output, retry on transient errors, have a fallback.
6. **Add caching.** Cache responses keyed on (model, prompt hash) to reduce
   cost and latency on repeated identical inputs.
7. **Evaluate.** Write a small eval suite — a set of input/expected-output
   pairs that runs in CI. Prompt changes are code changes; test them.

## Skeleton — Anthropic Claude (Python)

```python
# llm/client.py
import os
import anthropic
from pydantic import BaseModel

client = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"])

SYSTEM_PROMPT = """You are a precise data extractor. Given unstructured text,
extract the requested fields exactly as they appear. If a field is absent,
return null. Do not infer or hallucinate values."""


class ExtractedEntity(BaseModel):
    name: str | None
    date: str | None
    amount: float | None
    currency: str | None


def extract_entity(text: str, model: str = "claude-sonnet-4-5") -> ExtractedEntity:
    message = client.messages.create(
        model=model,
        max_tokens=512,
        system=SYSTEM_PROMPT,
        messages=[
            {
                "role": "user",
                "content": f"Extract the entity from this text:\n\n{text}\n\nRespond only with valid JSON matching this schema: {ExtractedEntity.model_json_schema()}"
            }
        ],
    )
    raw = message.content[0].text
    return ExtractedEntity.model_validate_json(raw)
```

## Skeleton — streaming response (FastAPI + SSE)

```python
# api/chat.py
from fastapi import APIRouter
from fastapi.responses import StreamingResponse
import anthropic
import json

router = APIRouter()
client = anthropic.Anthropic()


@router.post("/chat")
async def chat(prompt: str):
    def generate():
        with client.messages.stream(
            model="claude-sonnet-4-5",
            max_tokens=1024,
            messages=[{"role": "user", "content": prompt}],
        ) as stream:
            for text in stream.text_stream:
                yield f"data: {json.dumps({'text': text})}\n\n"
        yield "data: [DONE]\n\n"

    return StreamingResponse(generate(), media_type="text/event-stream")
```

## Skeleton — tool use / function calling

```python
# llm/tools.py
import anthropic
import json

client = anthropic.Anthropic()

tools = [
    {
        "name": "get_order_status",
        "description": "Look up the status of a customer order by order ID.",
        "input_schema": {
            "type": "object",
            "properties": {
                "order_id": {"type": "string", "description": "The order ID"}
            },
            "required": ["order_id"],
        },
    }
]


def run_with_tools(user_message: str) -> str:
    messages = [{"role": "user", "content": user_message}]

    while True:
        response = client.messages.create(
            model="claude-sonnet-4-5",
            max_tokens=1024,
            tools=tools,
            messages=messages,
        )

        if response.stop_reason == "end_turn":
            return response.content[0].text

        if response.stop_reason == "tool_use":
            tool_call = next(b for b in response.content if b.type == "tool_use")
            # Execute the tool
            result = dispatch_tool(tool_call.name, tool_call.input)
            # Feed result back
            messages += [
                {"role": "assistant", "content": response.content},
                {"role": "user", "content": [{"type": "tool_result", "tool_use_id": tool_call.id, "content": result}]},
            ]


def dispatch_tool(name: str, inputs: dict) -> str:
    if name == "get_order_status":
        return json.dumps({"order_id": inputs["order_id"], "status": "shipped"})
    raise ValueError(f"Unknown tool: {name}")
```

## Prompt engineering checklist

Before shipping any prompt:

- [ ] System prompt is a persona + constraints, not instructions (those go in user).
- [ ] Output format is explicitly specified with an example or JSON schema.
- [ ] Edge cases are listed: "if X is missing, return null" not "return your best guess".
- [ ] Tested with adversarial inputs (empty string, malformed data, very long text).
- [ ] Prompt is versioned — stored as a string constant, not scattered through business logic.
- [ ] Cost estimated: `input_tokens * price + output_tokens * price` at expected volume.

## Eval skeleton

```python
# tests/test_prompt.py
import pytest
from llm.client import extract_entity

CASES = [
    ("Invoice #1234 for $450.00 USD due 2025-06-01", {"amount": 450.0, "currency": "USD"}),
    ("No financial data here", {"amount": None, "currency": None}),
    ("€1,200 invoice", {"amount": 1200.0, "currency": "EUR"}),
]

@pytest.mark.parametrize("text,expected", CASES)
def test_extract_entity(text, expected):
    result = extract_entity(text)
    for key, val in expected.items():
        assert getattr(result, key) == val, f"Field {key}: expected {val}, got {getattr(result, key)}"
```

## Guardrails

- Follow [`../../shared/guidelines/security-baseline.md`](../../shared/guidelines/security-baseline.md):
  API keys from env only; never log user content that may be PII.
- Never send PII to an external model without user consent and a data-handling agreement.
- Set `max_tokens` explicitly — unbounded generation is a cost and latency risk.
- Validate structured outputs before using them downstream; models can return malformed JSON.
- Cache aggressively on the prompt hash to avoid redundant API calls.

## Anti-patterns

- Parsing JSON from raw model output with regex — use structured output APIs.
- Putting all prompt logic in a single mega-string — split system / user / few-shot examples.
- No fallback when the model fails — always handle `anthropic.APIError`, timeouts, and invalid output.
- Testing the prompt only manually — add even 5 representative eval cases to CI.
- Hardcoding the model string everywhere — define it as a config value.

## Extending

- **Multi-turn chat:** maintain a `messages: list[dict]` state and append each turn.
- **RAG:** embed user query → retrieve relevant chunks → inject as context in the user prompt.
- **Cost tracking:** wrap the client to log token counts per call; alert if daily cost spikes.
- **Model routing:** cheap model for classification, expensive model for generation, based on task complexity.
