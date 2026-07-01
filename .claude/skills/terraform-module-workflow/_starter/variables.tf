# ---------------------------------------------------------------------------
# Naming and tagging — required on every instance
# ---------------------------------------------------------------------------

variable "name" {
  description = "Unique name for this module instance. Used as a prefix for all resources."
  type        = string

  validation {
    condition     = can(regex("^[a-z0-9-]{3,40}$", var.name))
    error_message = "Name must be 3-40 chars, lowercase alphanumeric or hyphens."
  }
}

variable "environment" {
  description = "Deployment environment (e.g. dev, stage, prod)."
  type        = string

  validation {
    condition     = contains(["dev", "stage", "prod", "sandbox"], var.environment)
    error_message = "Environment must be one of: dev, stage, prod, sandbox."
  }
}

variable "tags" {
  description = "Additional tags / labels merged into the common tag map applied to every resource."
  type        = map(string)
  default     = {}
}

# ---------------------------------------------------------------------------
# Cloud-specific placement — pick the subset your module needs and delete the rest
# ---------------------------------------------------------------------------

variable "aws_region" {
  description = "AWS region. Leave null if this module doesn't target AWS."
  type        = string
  default     = null
}

variable "gcp_project" {
  description = "GCP project ID. Leave null if this module doesn't target GCP."
  type        = string
  default     = null
}

variable "gcp_region" {
  description = "GCP region."
  type        = string
  default     = null
}

variable "azure_subscription_id" {
  description = "Azure subscription GUID. Leave null if this module doesn't target Azure."
  type        = string
  default     = null
}

variable "azure_location" {
  description = "Azure location (e.g. eastus, westeurope)."
  type        = string
  default     = null
}

# ---------------------------------------------------------------------------
# Module-specific inputs — add yours below.
# ---------------------------------------------------------------------------

# variable "example_input" {
#   description = "Describe what it does and why it matters."
#   type        = string
# }
