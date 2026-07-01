---
name: database-migration
description: Use when the user asks to add, alter, or drop database tables or columns, create a migration script, or run schema changes safely against a live database.
allowed-tools: Read, Grep, Glob, Bash, Edit, Write
---

# database-migration

Safe, reviewable database schema migrations: write the migration, write
the rollback, verify idempotency, and walk through the deployment checklist.
Never runs `apply` or `upgrade` without the user's explicit instruction.

## When to invoke

- User asks to "add a column", "rename a table", "create an index", "drop a field".
- User asks to scaffold an Alembic or Django migration.
- User wants to backfill data as part of a schema change.
- User is debugging a failed migration or a drift between schema and ORM models.

## Method

1. **Read the ORM models and existing migrations.** Understand current schema,
   naming conventions, and the migration framework in use.
2. **Draft the forward migration** (schema change only; data backfill is separate).
3. **Draft the rollback** (`downgrade` in Alembic; `reverse` in Django). If the
   change is irreversible (drop column), say so explicitly.
4. **Check for data loss.** `DROP COLUMN` and `NOT NULL` additions without
   defaults are destructive. Present the user with the risk and a safe two-step
   alternative if needed.
5. **Verify idempotency.** Use `IF EXISTS` / `IF NOT EXISTS` guards.
6. **Walk through the deployment checklist** (see below).
7. **Stop.** Do not run `alembic upgrade head` or `python manage.py migrate`
   unless the user explicitly asks.

## Skeletons

### Alembic (SQLAlchemy)

```python
"""add status column to orders

Revision ID: abc123def456
Revises: 789ghi012jkl
Create Date: 2025-01-01 00:00:00.000000
"""
from alembic import op
import sqlalchemy as sa

revision = "abc123def456"
down_revision = "789ghi012jkl"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column(
        "orders",
        sa.Column(
            "status",
            sa.String(length=32),
            nullable=False,
            server_default="pending",  # required for NOT NULL on existing rows
        ),
    )
    op.create_index("ix_orders_status", "orders", ["status"])


def downgrade() -> None:
    op.drop_index("ix_orders_status", table_name="orders")
    op.drop_column("orders", "status")
```

### Django

```python
# migrations/0042_add_status_to_order.py
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [("myapp", "0041_previous_migration")]

    operations = [
        migrations.AddField(
            model_name="order",
            name="status",
            field=models.CharField(max_length=32, default="pending"),
            preserve_default=False,
        ),
    ]
```

### Raw SQL (with idempotency guard)

```sql
-- up.sql
ALTER TABLE orders
    ADD COLUMN IF NOT EXISTS status VARCHAR(32) NOT NULL DEFAULT 'pending';

CREATE INDEX IF NOT EXISTS ix_orders_status ON orders (status);

-- down.sql
DROP INDEX IF EXISTS ix_orders_status;
ALTER TABLE orders DROP COLUMN IF EXISTS status;
```

## Safe patterns for destructive changes

### Adding NOT NULL without downtime (three-step)

```
Step 1: Add nullable column with no default → deploy
Step 2: Backfill existing rows → run as a job
Step 3: Add NOT NULL constraint → deploy
```

### Renaming a column without downtime

```
Step 1: Add new column → deploy
Step 2: Dual-write (write both old + new) → deploy
Step 3: Backfill old → null out old → deploy
Step 4: Drop old column → deploy
```

## Deployment checklist

- [ ] Migration is reviewed and approved.
- [ ] Rollback script tested on a copy of production data.
- [ ] Large tables: estimate lock time. Index `CONCURRENTLY` on Postgres.
- [ ] Alembic: `alembic current` shows expected revision before upgrade.
- [ ] Django: `./manage.py showmigrations` — all preceding are applied.
- [ ] Staged rollout: apply to staging first; compare row counts before/after.
- [ ] On-call team notified if change touches tables with > 1M rows or active writes.

## Guardrails

- Follow [`../../shared/guidelines/security-baseline.md`](../../shared/guidelines/security-baseline.md):
  never include credentials in migration files.
- Never `DROP TABLE` or `DROP COLUMN` without a confirmed backup and explicit user instruction.
- Don't add a `NOT NULL` column to a table with existing rows without a `DEFAULT` or a backfill plan.
- Migration files are append-only once merged — never edit a migration that has run in production.

## Anti-patterns

- Running migrations inside a deployment that also restarts the application — do migrations first, then deploy.
- One giant migration that adds columns, backfills data, and drops old columns in a single transaction.
- Skipping `downgrade()` — every migration needs a tested rollback.
- `op.execute("ALTER TABLE ...")` without `IF EXISTS` guards in shared environments.
