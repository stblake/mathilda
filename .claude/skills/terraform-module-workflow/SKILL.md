---
name: terraform-module-workflow
description: Use when the user asks to create, modify, or refactor Terraform — a new module, adding a resource to an existing module, or splitting a monolithic config into reusable modules. Produces portable, multi-cloud-aware infrastructure.
---

# terraform-module-workflow

The default workflow for authoring or evolving Terraform. Opinionated toward
portable module shapes so projects can retarget clouds without a rewrite.

**Ships a runnable starter.** [`_starter/`](_starter/) is a ready-to-copy
module directory with the canonical layout, tagging local, variable
validation, and a working `examples/basic/` consumer. Start there:

```bash
cp -r skills/terraform-module-workflow/_starter my-project/modules/<name>
cd my-project/modules/<name>/examples/basic
cp terraform.tfvars.example terraform.tfvars
terraform init && terraform plan
```

## When to invoke

- User asks to "create a Terraform module for X" or "scaffold infra for X".
- User points at an empty `modules/<name>/` directory.
- User wants to refactor a sprawling `main.tf` into reusable modules.
- User asks to add a resource (bucket, queue, DB, function) to existing IaC.

## Method

1. **Read the project.** Are there existing modules? What provider versions?
   What backend is configured? Match conventions rather than impose new ones.
2. **Confirm the cloud target.** AWS, GCP, or Azure. A module targets one
   cloud; portability comes from the module *interface*, not from a switch.
3. **Scaffold to the canonical layout** (below) unless the project differs.
4. **Fill the interface.** Inputs, outputs, tags, naming — uniform across
   all our modules.
5. **Add a minimal example** under `examples/basic/` so a consumer can copy
   it and get a working apply.
6. **Validate.** `terraform fmt` and `terraform validate` on the example.
   Don't `terraform apply` from the skill; the user runs apply.

## Canonical module layout

```
modules/<name>/
  README.md           # purpose, inputs, outputs, one example
  main.tf             # resources
  variables.tf        # input contract
  outputs.tf          # output contract
  versions.tf         # terraform + provider pins
  providers.tf        # provider config (region/project/subscription as vars)
  examples/
    basic/
      main.tf
      terraform.tfvars.example
```

## Conventions to enforce

See [`../../shared/guidelines/multi-cloud.md`](../../shared/guidelines/multi-cloud.md) and
[`../../shared/guidelines/security-baseline.md`](../../shared/guidelines/security-baseline.md).

**`versions.tf`** — pin `terraform` to a compatible range; pin each provider
with `~>` on the minor.

**`variables.tf`** — every variable has `description` and `type`. Required
inputs have no default. Validation blocks on enumerations and formats.

**`outputs.tf`** — IDs, ARNs, self-links. Mark `sensitive = true` where
appropriate.

**`main.tf`** — no literal region / project / subscription; always a
variable. Tags / labels merged from a common local map. Encryption at rest
on by default for stateful resources.

**`providers.tf`** — no credentials. Alias when the module uses multiple
regions.

## Skeletons

```hcl
# versions.tf (AWS example)
terraform {
  required_version = ">= 1.6, < 2.0"
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }
}
```

```hcl
# variables.tf
variable "name" {
  description = "Unique name for this module instance. Used as a prefix for all resources."
  type        = string
  validation {
    condition     = can(regex("^[a-z0-9-]{3,40}$", var.name))
    error_message = "Name must be 3-40 chars, lowercase alphanumeric or hyphens."
  }
}

variable "tags" {
  description = "Tags applied to every resource created by this module."
  type        = map(string)
  default     = {}
}
```

```hcl
# main.tf — tagging local
locals {
  common_tags = merge(var.tags, {
    "managed-by" = "terraform"
    "module"     = "<name>"
  })
}
```

## Guardrails

- **No backend config inside modules.** Backends belong to root modules.
- **No `terraform apply` from the skill.** Scaffolding and validation only.
- **No real account IDs, project IDs, subscription GUIDs** in examples —
  use `REPLACE_ME`.
- **No client names** anywhere.

## Extending

- Cloud pack: mirror the module under each cloud when a team truly needs
  multi-cloud parity (e.g. `modules/object-store/aws/`, `/gcp/`, `/azure/`).
- Compliance pack: overlay modules add audit logging, encryption keys, VPC
  isolation for regulated environments.
