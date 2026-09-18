args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  crops <- args
} else {
  crops <- c("Lettuce", "Strawberry", "Spinach", "Asparagus", "Cabbage")
}
