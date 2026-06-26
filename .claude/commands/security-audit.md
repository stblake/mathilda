---
description: Run a security audit on the current branch — surfaces credential leaks, injection risks, overly broad permissions, and dependency vulnerabilities.
argument-hint: [--deep]
allowed-tools: Read, Grep, Glob, Bash
---

Run a security audit on the changes in the current branch and the overall
codebase. Pass `--deep` as `$ARGUMENTS` to scan the full repo, not just the diff.

## Step 1 — Scope

```bash
# Changes on this branch vs main
git diff main...HEAD --name-only
```

If `--deep` was passed, audit the full repo. Otherwise, limit to changed files.

## Step 2 — Credential and secret scan

```bash
# Use detect-secrets if available
detect-secrets scan --all-files --baseline .secrets.baseline 2>/dev/null || \
  grep -r -n \
    -e "api_key\s*=\s*['\"]" \
    -e "secret\s*=\s*['\"]" \
    -e "password\s*=\s*['\"]" \
    -e "token\s*=\s*['\"]" \
    -e "AKIA[0-9A-Z]{16}" \
    --include="*.py" --include="*.ts" --include="*.js" --include="*.env" \
    --include="*.yaml" --include="*.json" \
    . | grep -v ".git" | grep -v "test_" | grep -v "__tests__"
```

Flag any hardcoded credentials. Check `.env` files are in `.gitignore`.

## Step 3 — Input validation audit

In changed files, check:

```bash
# Python: look for raw string interpolation into SQL
grep -n "f\".*SELECT\|f\".*INSERT\|f\".*UPDATE\|f\".*DELETE\|% .*cursor\|format.*cursor" \
  $(git diff main...HEAD --name-only | grep "\.py$")

# Check for shell injection
grep -n "subprocess.*shell=True\|os\.system(" \
  $(git diff main...HEAD --name-only | grep "\.py$")

# JS/TS: look for XSS risks
grep -n "dangerouslySetInnerHTML\|innerHTML\s*=" \
  $(git diff main...HEAD --name-only | grep -E "\.(tsx?|jsx?)$")
```

## Step 4 — Dependency vulnerability scan

```bash
# Python
pip audit 2>/dev/null || safety check 2>/dev/null || \
  echo "Install pip-audit: pip install pip-audit"

# Node
npm audit --audit-level=high 2>/dev/null || \
  echo "Run: npm audit"
```

List any HIGH or CRITICAL CVEs found.

## Step 5 — IAM and network exposure (Terraform)

If Terraform files changed:

```bash
# Overly broad IAM
grep -rn '"Action": "\*"\|"Resource": "\*"' \
  $(git diff main...HEAD --name-only | grep "\.tf$")

# Open security groups
grep -rn 'cidr_blocks.*0\.0\.0\.0/0\|cidr.*= \["0\.0\.0\.0/0"\]' \
  $(git diff main...HEAD --name-only | grep "\.tf$")

# Public S3 buckets
grep -rn 'acl.*=.*"public\|block_public_acls.*=.*false' \
  $(git diff main...HEAD --name-only | grep "\.tf$")
```

## Step 6 — Report

```
## Security Audit — <branch>

### Blockers (must fix before merge)
- [credential|injection|permission|cve] <file:line> — <description> → <fix>

### Warnings (should fix)
- [type] <file:line> — <description> → <suggestion>

### Informational
- <observation with no immediate action>

### Scan coverage
- Files audited: N
- Tools run: detect-secrets / pip-audit / npm audit / manual grep
- Limitations: <anything not covered>
```

## References

- [`../shared/guidelines/security-baseline.md`](../shared/guidelines/security-baseline.md)
- Command: [`review-pr`](review-pr.md) — includes security as one axis.
