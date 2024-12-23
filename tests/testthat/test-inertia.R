# ---- DATA PREPARATION ----
data("iris")
df0 <- select(iris, -Species)
df <- scale(df0)

# Perform hierarchical clustering using Ward's method on Euclidean distance
clustering_ward <- dist(df) %>% getCAH(3, .)
cls <- getClusterPerPart(c = clustering_ward)  # Extract cluster partitions
cl_full <- getClusterPerPart(c = clustering_ward, n = 150)  # Partition with 150 clusters

# ---- TEST: RELATIVE BETWEEN-CLUSTER VARIANCE ----
test_that("betweenPerPart", {
  res <- getRelativeBetweenPerPart(d = df, cl = cls)

  # Expected relative variances
  round(res, 5) %>%
    expect_identical(c(61.87981, 75.18729, 80.52696, 84.16555, 85.66615))
  # getRelativeBetweenPerPart(d = df0, cl = cls) %>%
  #   expect_equal(res)

  # Incremental differences in variance
  getBetweenDifferences(res) %>%
    round(6) %>%
    expect_identical(c(61.879806, 13.307485, 5.339665, 3.638599, 1.500597))

  # Full partition should explain 100% variance
  expect_identical(getRelativeBetweenPerPart(df, cl_full, n = 150)[149], 100)
})

# ---- TEST: WITHIN-CLUSTER VARIANCE ----
test_that("withinPerCluster", {
  # Calculate within-cluster variance using the original unscaled data
  res <- getRelativeWithinPerCluster(cls = cls, d = df0)

  # Verify consistency between original and scaled data
  getRelativeWithinPerCluster(cls = cls, d = df) %>%
    expect_equal(res)

  # Confirm within-cluster variance sums to 1
  apply(res, 1, function(x) sum(x, na.rm = TRUE)) %>%
    round(6) %>%
    unique() %>%
    expect_equal(1)

  # Dimension check: 5 partitions, max 6 clusters
  expect_equal(dim(res), c(5, 6))

  # Specific value check for the first partition
  res[1, ] %>%
    unique() %>%
    round(7) %>%
    expect_identical(c(0.6733333, 0.3266667, NA))
  # TODO: after getRelativeWithinPerCluster optimization, test for n = 149
})

# ---- TEST: DISTANCE PER VARIABLE ----
test_that("distPerVariable", {
  # Calculate distances per variable for the first partition
  res <- getDistPerVariable(d = df0, cl = cls[[1]])

  # Distances should be symmetric and sum to zero
  expect_equal(res[1, ], -res[2, ])
  apply(res, 2, sum) %>%
    round(7) %>%
    as.vector() %>%
    expect_equal(rep(0, 4))

  # Verify the matrix has the expected number of rows (2 clusters)
  expect_equal(nrow(res), 2)

  # Test full partition (150 clusters)
  res <- getDistPerVariable(d = df0, cl = cl_full[[149]])
  apply(res, 2, sum) %>%
    round(7) %>%
    as.vector() %>%
    expect_equal(rep(0, 4))
  expect_equal(nrow(res), 150)

  # Ensure consistency with scaled data
  getDistPerVariable(d = df, cl = cl_full[[149]]) %>%
    expect_equal(res)
})

# ---- TEST: VARIABLE CONTRIBUTIONS ----
test_that("getCtrVar", {
  # Calculate variable contributions to clustering for 3 clusters (original and scaled data)
  res <- getCtrVar2(cl = cls[[2]], d = df0, k = 3, scale = FALSE)
  res2 <- getCtrVar2(cl = cls[[2]], d = df, k = 3)

  # Ensure results differ between scaled and unscaled data
  unlist(res2) %>%
    round(7) %>%
    as.vector() %>%
    setdiff(unlist(res) %>% as.numeric() %>% round(7)) %>%
    length() %>%
    expect_gt(0)

  # Contributions should sum to 1
  apply(res2, 2, sum) %>%
    round(6) %>%
    as.vector() %>%
    expect_equal(rep(1, 4))

  # Check consistency between two contribution calculation methods
  res2 <- getCtrVar(cl = cls[[2]], d = df0, k = 3)
  expect_equal(res, res2)
  getCtrVar(cl = cls[[2]], d = df, k = 3) %>%
    expect_equal(res2)

  # Check specific contribution values for the first variable
  res[, 1] %>%
    round(6) %>%
    as.vector() %>%
    expect_equal(c(0.32801820, 0.02882824, 0.34358365) %>% round(6))

  # Validate results for full partition (150 clusters)
  getCtrVar2(cl = cl_full[[149]], d = df0, k = 150, scale = FALSE) %>%
    dim() %>%
    expect_equal(c(150, 4))
})

# ---- TEST: PARTITION-WISE DISTANCE CONTRIBUTIONS ----
test_that("getPdis", {
  # Calculate partition distance for 3 clusters
  res <- getPdis(cl = cls[[2]], d = df, k = 3)

  # Partition distance should match sum of variable contributions
  getCtrVar(cl = cls[[2]], d = df, k = 3) %>%
    apply(2, sum) %>%
    expect_equal(res)

  # Test partition distances for multiple partitions
  res2 <- getPdisPerPartition(cls = cls, d = df)
  expect_equal(res2[2, ], res)

  # With n = 7
  res2 <- getPdisPerPartition(cls = cl_full, d = df, n = 7)
  expect_equal(dim(res2), c(6, 4))

  # Specific value checks
  res2[, 1] %>%
    round(7) %>%
    as.vector() %>%
    expect_equal(c(0.4871557, 0.7004301, 0.7901667, 0.8192104, 0.8425190, 0.8660134))
})
