data("iris")
df0 <- select(iris, -Species)
df <- scale(df0)

# Perform hierarchical clustering using Ward's method on Euclidean distance
dist <- dist(df)
clustering_ward <- getCAH(3, dist)

# Extract cluster partitions (up to 10 clusters)
cls <- getClusterPerPart(c = clustering_ward, n = 10)

# Compute the gap statistic for cluster evaluation (B = 10 resamples)
gap_res <- getGapPerPart(d = df0, c = clustering_ward, B = 10, n = 10)

# ---- TEST GAP STATISTIC OUTPUT ----
test_that("gapPerPart", {
  # Ensure gap statistic result is of the correct class
  expect_s3_class(gap_res, "clusGap")

  # Verify that the result contains 10 rows (one per cluster partition)
  expect_equal(nrow(gap_res$Tab), 10)
})

# ---- TEST OPTIMAL NUMBER OF CLUSTERS ----
test_that("gapBest", {
  res <- getGapBest(gap_res)

  # Optimal number of clusters should be between 1 and 10
  expect_gte(res, 1)
  expect_lte(res, 10)
})

# ---- TEST CLUSTER SUMMARY OUTPUT ----
test_that("printSummary", {
  # Calculate clustering metrics
  between <- getRelativeBetweenPerPart(d = df, cl = cls, n = 10)  # Between-cluster variance
  diff <- getBetweenDifferences(between)  # Incremental variance differences
  silhouette <- getSilhouettePerPart(cls, dist) %>%
    getMeanSilhouettePerPart()  # Average silhouette width

  # Generate summary table of clustering performance
  res <- printSummary(between, diff, silhouette)

  # Ensure the summary table has 9 rows and 4 columns
  expect_equal(dim(res), c(9, 4))

  # Validate the last row of the summary
  as.vector(res[9, ]) %>%
    round(7) %>%
    expect_equal(c(91.119214,  1.288824,  8.880786,  0.316112))

  # Recompute summary including gap statistics
  res <- printSummary(between, diff, silhouette, gap_res)

  # Ensure the summary now has 6 columns (additional gap statistics)
  expect_equal(dim(res), c(9, 6))

  # Check that between-cluster variance, differences, and silhouette match
  expect_equal(round(as.vector(res[, 1]), 7), between)
  expect_equal(round(as.vector(res[, 2]), 7), diff)
  expect_equal(round(as.vector(res[, 4]), 7), round(silhouette, 7))
})

