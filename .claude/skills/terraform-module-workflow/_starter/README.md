# `<module-name>`

> Replace this title and description after copying the starter.

One-paragraph description: what the module provisions, which cloud it
targets, what consumers should expect.

## Usage

```hcl
module "<module-name>" {
  source = "../../modules/<module-name>"

  name        = "my-instance"
  environment = "dev"
  tags = {
    Project = "example"
    Owner   = "team@example.com"
  }

  # cloud-specific
  aws_region = "us-east-1"

  # module-specific inputs
  # ...
}
```

A runnable example lives under [`examples/basic/`](examples/basic/).

## Inputs

See [`variables.tf`](variables.tf). Required:

- `name` — resource name prefix.
- `environment` — one of `dev`, `stage`, `prod`, `sandbox`.

## Outputs

See [`outputs.tf`](outputs.tf).

## Neutral alternatives on other clouds

| This module's resource | AWS | GCP | Azure |
| ---------------------- | --- | --- | ----- |
| (describe)             |     |     |       |

See [`shared/guidelines/multi-cloud.md`](../../../shared/guidelines/multi-cloud.md).

## Conventions enforced

- Every resource tagged via `local.common_tags`.
- Encryption at rest on stateful resources.
- No `0.0.0.0/0` on non-public ports.
- Variables for region / project / subscription — no literals.

## Testing

```bash
cd examples/basic
cp terraform.tfvars.example terraform.tfvars
# edit terraform.tfvars
terraform init
terraform plan
```

Do not `terraform apply` from a clean clone without reviewing the plan.
