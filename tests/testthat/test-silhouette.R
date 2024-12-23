# ---- DATA PREPARATION ----
data("iris")
df0 <- select(iris, -Species)
df <- scale(df0)

# Perform hierarchical clustering using Ward's method on Euclidean distance
dist <- dist(df)
clustering_ward <- getCAH(3, dist)

# Extract cluster partitions
cls <- getClusterPerPart(c = clustering_ward)  # Default partition
cl_full <- getClusterPerPart(c = clustering_ward, n = 150)  # Full partitioning (150 clusters)

# Calculate silhouette metrics for partitions
silhouette_k <- getSilhouette(cls[[1]], dist)  # Silhouette for the first partition
silhouette_148 <- getSilhouette(cl_full[[148]], dist)  # Silhouette for 148 clusters
silhouettes <- getSilhouettePerPart(cl_full[-149], dist)  # Silhouette for all but the last partition

# ---- SILHOUETTE TEST FUNCTION ----
test_silhouette <- function(x) {
  # Check if result is a silhouette object
  expect_s3_class(x, "silhouette")

  # Ensure silhouette values are within valid range [-1, 1]
  expect_gte(min(x[, 3]), -1)
  expect_lte(max(x[, 3]), 1)

  # Validate that silhouette length corresponds to 150 data points
  expect_length(x[, 1], 150)
}

# ---- TEST INDIVIDUAL SILHOUETTE CALCULATION ----
test_that("silhouette", {
  # Ensure clustering result only contains 2 clusters (1 and 2)
  expect_equal(unique(silhouette_k[, 1]), seq(2))

  # Validate first six silhouette values for partition 1
  head(silhouette_k[, 3]) %>%
    as.numeric() %>%
    round(7) %>%
    expect_equal(c(0.7746877, 0.7744370, 0.7739863, 0.7699909, 0.7691272, 0.7669716))

  # Perform general silhouette checks
  test_silhouette(silhouette_k)
  test_silhouette(silhouette_148)
})

# ---- TEST SILHOUETTE CONSISTENCY ACROSS PARTITIONS ----
test_that("silhouettePerPart", {
  # Verify partition 1 in silhouettes matches silhouette_k
  expect_equal(silhouettes[[1]], silhouette_k)

  # Validate partition 148 matches the 148-cluster silhouette
  expect_equal(silhouettes[[148]], silhouette_148)
})

# ---- TEST MEAN SILHOUETTE VALUES ----
test_that("meanSilhouette", {
  # Calculate mean silhouette for all partitions
  # Ensure mean silhouette for partition 148 matches the expected value
  getMeanSilhouettePerPart(silhouettes) %>%
    .[148] %>%
    expect_equal(silhouette_148[, 3] %>% mean())
})
