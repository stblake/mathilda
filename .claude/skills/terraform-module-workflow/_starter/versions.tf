terraform {
  required_version = ">= 1.6, < 2.0"

  # Keep only the providers the module actually uses.
  # Delete the blocks for clouds this module doesn't target.
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
    google = {
      source  = "hashicorp/google"
      version = "~> 5.0"
    }
    azurerm = {
      source  = "hashicorp/azurerm"
      version = "~> 3.0"
    }
  }
}
