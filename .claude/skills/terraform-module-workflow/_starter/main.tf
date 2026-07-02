# ---------------------------------------------------------------------------
# Common tags applied to every resource in this module.
# Consumers get their overlay via var.tags; module identity stays fixed.
# ---------------------------------------------------------------------------

locals {
  module_name = "REPLACE_ME_MODULE_NAME"

  common_tags = merge(
    var.tags,
    {
      Name        = var.name
      Environment = var.environment
      ManagedBy   = "terraform"
      Module      = local.module_name
    },
  )
}

# ---------------------------------------------------------------------------
# Resources go below. Prefer editing and copying existing resource blocks
# over introducing new patterns. Follow the security baseline:
#
#   - private by default for storage
#   - encryption at rest on stateful resources
#   - no 0.0.0.0/0 on non-public ports
#   - least-privilege IAM
#
# See: shared/guidelines/security-baseline.md and multi-cloud.md.
# ---------------------------------------------------------------------------

# Example (delete and replace):
#
# resource "aws_s3_bucket" "this" {
#   bucket = "${var.name}-${var.environment}"
#   tags   = local.common_tags
# }
#
# resource "aws_s3_bucket_server_side_encryption_configuration" "this" {
#   bucket = aws_s3_bucket.this.id
#   rule {
#     apply_server_side_encryption_by_default {
#       sse_algorithm = "AES256"
#     }
#   }
# }
#
# resource "aws_s3_bucket_public_access_block" "this" {
#   bucket                  = aws_s3_bucket.this.id
#   block_public_acls       = true
#   block_public_policy     = true
#   ignore_public_acls      = true
#   restrict_public_buckets = true
# }
