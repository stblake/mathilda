---
name: monitoring-observability
description: Use when the user asks to add logging, metrics, tracing, health checks, or alerting to a service — or when a service lacks visibility into what it's doing in production.
allowed-tools: Read, Grep, Glob, Bash, Edit, Write
---

# monitoring-observability

Add production-grade observability to a service: structured logging, request
metrics, distributed tracing, health/readiness endpoints, and a basic alerting
contract. Ships wiring code; choose the backend that matches the project.

## When to invoke

- User says "add logging", "set up metrics", "instrument this service".
- Service has `print()` statements instead of structured logs.
- There is no `/health` or `/ready` endpoint.
- User is debugging a production issue and has no traces or metrics to look at.
- User is adding a new service and wants observability from day one.

## Method

1. **Read the project.** What framework? What observability stack is already
   in use (Datadog, OpenTelemetry, CloudWatch, Prometheus)?
2. **Add structured logging.** Replace `print` with `structlog` or `logging`
   with a JSON formatter. Add request ID to every log line.
3. **Add metrics.** Request count, latency histogram, error rate — the RED
   method (Rate, Errors, Duration). Use `prometheus_client` or the project's
   existing SDK.
4. **Add tracing.** Instrument the entry point and outgoing HTTP/DB calls with
   OpenTelemetry spans. Propagate context across service boundaries.
5. **Add health endpoints.** `/health` (liveness) and `/ready` (readiness,
   checks DB + dependencies).
6. **Document alert thresholds.** At minimum: error rate > 1%, p99 latency > SLO.

## Logging skeleton (structlog + FastAPI)

```python
# app/logging_config.py
import logging
import structlog

def configure_logging(log_level: str = "INFO") -> None:
    structlog.configure(
        processors=[
            structlog.contextvars.merge_contextvars,
            structlog.processors.add_log_level,
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.processors.JSONRenderer(),
        ],
        wrapper_class=structlog.make_filtering_bound_logger(
            getattr(logging, log_level.upper())
        ),
        context_class=dict,
        logger_factory=structlog.PrintLoggerFactory(),
    )

log = structlog.get_logger()
```

```python
# app/middleware.py
import uuid
import time
import structlog
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request

log = structlog.get_logger()


class RequestLoggingMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next):
        request_id = str(uuid.uuid4())
        structlog.contextvars.bind_contextvars(request_id=request_id)
        start = time.perf_counter()
        response = await call_next(request)
        duration_ms = (time.perf_counter() - start) * 1000
        log.info(
            "request",
            method=request.method,
            path=request.url.path,
            status=response.status_code,
            duration_ms=round(duration_ms, 2),
        )
        structlog.contextvars.clear_contextvars()
        response.headers["X-Request-ID"] = request_id
        return response
```

## Metrics skeleton (Prometheus)

```python
# app/metrics.py
from prometheus_client import Counter, Histogram, generate_latest, CONTENT_TYPE_LATEST
from fastapi import APIRouter
from fastapi.responses import Response

REQUEST_COUNT = Counter(
    "http_requests_total",
    "Total HTTP requests",
    ["method", "path", "status"],
)
REQUEST_DURATION = Histogram(
    "http_request_duration_seconds",
    "HTTP request latency",
    ["method", "path"],
    buckets=[0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0, 2.5],
)

router = APIRouter()

@router.get("/metrics", include_in_schema=False)
def metrics():
    return Response(generate_latest(), media_type=CONTENT_TYPE_LATEST)
```

## Health endpoints skeleton

```python
# app/health.py
from fastapi import APIRouter
from pydantic import BaseModel
from app.db import get_db_status  # implement per project

router = APIRouter()

class HealthResponse(BaseModel):
    status: str
    checks: dict[str, str]

@router.get("/health", response_model=HealthResponse)
async def liveness():
    return {"status": "ok", "checks": {}}

@router.get("/ready", response_model=HealthResponse)
async def readiness():
    checks = {"db": await get_db_status()}
    status = "ok" if all(v == "ok" for v in checks.values()) else "degraded"
    return {"status": status, "checks": checks}
```

## OpenTelemetry tracing skeleton

```python
# app/tracing.py
from opentelemetry import trace
from opentelemetry.sdk.trace import TracerProvider
from opentelemetry.sdk.trace.export import BatchSpanProcessor
from opentelemetry.exporter.otlp.proto.grpc.trace_exporter import OTLPSpanExporter

def configure_tracing(service_name: str, otlp_endpoint: str) -> None:
    provider = TracerProvider()
    exporter = OTLPSpanExporter(endpoint=otlp_endpoint)
    provider.add_span_processor(BatchSpanProcessor(exporter))
    trace.set_tracer_provider(provider)

tracer = trace.get_tracer(__name__)

# Usage in business logic:
# with tracer.start_as_current_span("fetch-user") as span:
#     span.set_attribute("user.id", user_id)
#     user = await db.get_user(user_id)
```

## Alert thresholds to document

```yaml
# alerts.yaml (Prometheus AlertManager format — adapt to your backend)
groups:
  - name: service-slos
    rules:
      - alert: HighErrorRate
        expr: rate(http_requests_total{status=~"5.."}[5m]) / rate(http_requests_total[5m]) > 0.01
        for: 2m
        labels:
          severity: critical
        annotations:
          summary: "Error rate > 1% for 2 minutes"

      - alert: HighLatency
        expr: histogram_quantile(0.99, rate(http_request_duration_seconds_bucket[5m])) > 1.0
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "p99 latency > 1s for 5 minutes"
```

## Guardrails

- Follow [`../../shared/guidelines/security-baseline.md`](../../shared/guidelines/security-baseline.md):
  never log PII, credentials, or raw request bodies.
- Structured logs only — no `print()`, no unstructured `logging.info("thing: " + str(x))`.
- `/health` must respond < 100ms. Don't put slow checks there — use `/ready`.
- Expose `/metrics` only on an internal port (not public-facing).
- Never log the `Authorization` header or cookie values.

## Anti-patterns

- Logging at `DEBUG` level in production by default — sets log level from config.
- One log line per SQL query without sampling — use `SQLALCHEMY_ECHO=false` in prod.
- Health checks that always return 200 regardless of DB state.
- Tracing spans that capture full request/response bodies — size and PII risk.
- Metrics with unbounded cardinality labels (e.g., `user_id` as a label).

## Extending

- **Cloud-native:** replace Prometheus with CloudWatch Metrics (AWS),
  Cloud Monitoring (GCP), or Azure Monitor.
- **Log aggregation:** ship structured JSON to Datadog, Loki, or CloudWatch Logs
  via a sidecar or OTLP exporter.
- **SLO tracking:** wrap `REQUEST_DURATION` with an SLO burn-rate alert.
- **Error tracking:** integrate `sentry-sdk` alongside structured logging.
