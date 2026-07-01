# Provider configuration. Keep the block(s) for clouds this module uses
# and delete the rest.
#
# NEVER put credentials here. Rely on:
#   - AWS: AWS_ACCESS_KEY_ID / AWS_PROFILE / instance profiles
#   - GCP: GOOGLE_APPLICATION_CREDENTIALS or workload identity
#   - Azure: az login / AZURE_CLIENT_ID / workload identity
#
# Alias providers when the module spans multiple regions.

# provider "aws" {
#   region = var.aws_region
# }
#
# provider "google" {
#   project = var.gcp_project
#   region  = var.gcp_region
# }
#
# provider "azurerm" {
#   subscription_id = var.azure_subscription_id
#   features {}
# }
