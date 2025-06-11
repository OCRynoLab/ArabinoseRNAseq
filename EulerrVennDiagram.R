# Install the Euler package if you don't already have it
if (!requireNamespace("eulerr", quietly = TRUE)) {
  install.packages("eulerr")

# Load the Eulerr package
library(eulerr)

# Define the sizes of sets and their overlaps
data <- c(
  Set1 = 573,          # Unique elements in Set1
  Set2 = 276,           # Unique elements in Set2
  "Set1&Set2" = 467,    # Overlap between Set1 and Set2
  Set3 = 76,           # Unique elements in Set3
  Set4 = 482,
  "Set1&Set3" = 196,    # Overlap between Set1 and Set3
  "Set2&Set3" = 103,    # Overlap between Set2 and Set3
  "Set1&Set4" = 542,
  "Set3&Set4" = 200,
  "Set2&Set4" = 246,
  "Set1&Set2&Set3" = 55, # Overlap between all three sets
  "Set2&Set3&Set4" = 45, # Overlap between all three sets
  "Set1&Set3&Set4" = 103, # Overlap between all three sets
  "Set1&Set2&Set4" = 134, # Overlap between all three sets
  "Set1&Set2&Set3&Set4" = 28 # Overlap between all three sets
)
fit <- euler(data)
plot(fit, 
     fills = c("green", "blue", "red","yellow"),  # Colors for the sets
     labels = TRUE,                      # Show set labels
     shape = "ellipse",
     main = "Scaled Venn Diagram")