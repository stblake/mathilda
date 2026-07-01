# Minimal example of consuming the module.
# Copy terraform.tfvars.example to terraform.tfvars, fill it in, then:
#   terraform init
#   terraform plan

module "example" {
  source = "../../"

  name        = var.name
  environment = var.environment

  aws_region = var.aws_region

  tags = {
    Project = "starter-example"
  }
}

variable "name" {
  type    = string
  default = "starter-example"
}

variable "environment" {
  type    = string
  default = "dev"
}

variable "aws_region" {
  type    = string
  default = "us-east-1"
}

# output "bucket_id" {
#   value = module.example.bucket_id
# }
