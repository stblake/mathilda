# Export resource IDs / ARNs / self-links that consumers will need to
# reference from downstream modules. Mark sensitive = true on anything
# that shouldn't appear in console output or state dumps.

# Example (delete and replace):
#
# output "bucket_id" {
#   description = "ID of the created object storage bucket."
#   value       = aws_s3_bucket.this.id
# }
#
# output "bucket_arn" {
#   description = "ARN of the created object storage bucket."
#   value       = aws_s3_bucket.this.arn
# }
#
# output "connection_string" {
#   description = "Database connection string. Sensitive."
#   value       = "postgresql://${local.username}:${local.password}@${aws_db_instance.this.endpoint}/${local.db_name}"
#   sensitive   = true
# }
