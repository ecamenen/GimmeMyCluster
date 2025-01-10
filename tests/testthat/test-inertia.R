# ---- DATA PREPARATION ----
data("iris")
df0 <- select(iris, -Species)
df <- scale(df0)

# Perform hierarchical clustering using Ward's method on Euclidean distance
clustering_ward <- dist(df) %>% run_ahc()
cls <- getClusterPerPart(clustering_ward)  # Extract cluster partitions
cl_full <- getClusterPerPart(clustering_ward, max_cluster = 150)  # Partition with 150 clusters

clustering_kmeans <- run_clustering(data = df, method = 2)
cl_kmeans <- extract_clusters(clustering_kmeans)$clusters

# ---- TEST: RELATIVE BETWEEN-CLUSTER VARIANCE ----
test_that("betweenPerPart", {
  res <- calculate_inertia(df0, cl = cls)

  # Expected relative variances
  round(res, 5) %>%
    expect_identical(c(74.82947, 85.42829, 89.79747, 90.72825, 92.20517))
  # calculate_inertia(d = df0, cl = cls) %>%
  #   expect_equal(res)

  # Incremental differences in variance
  calculate_between_diff(res) %>%
    round(7) %>%
    expect_identical(c(74.8294683, 10.5988240,  4.3691747,  0.9307804,  1.4769251))

  # Full partition should explain 100% variance
  expect_identical(calculate_inertia(df, cl = cl_full, max_cluster = 150)[149], 100)

  res <- calculate_between_inertia(df, cl = cl_kmeans)
  expect_equal(res[1], clustering_kmeans$result[[1]]$betweenss)
})

# ---- TEST: WITHIN-CLUSTER VARIANCE ----
test_that("withinPerCluster", {
  # Calculate within-cluster variance using the original unscaled data
  res <- calculate_within_inertia(df0, cl = cls)

  # Verify consistency between original and scaled data
  # calculate_within_inertia(df, cl = cls) %>%
  #   expect_equal(res)

  # Confirm within-cluster variance sums to 1
  # select(res, -1) %>%
  #   apply(1, function(x) sum(x, na.rm = TRUE)) %>%
  #   round(6) %>%
  #   unique() %>%
  #   expect_equal(100)

  # Dimension check: 5 partitions, max 6 clusters + rownames
  expect_equal(dim(res), c(5, 7))

  # Specific value check for the first partition
  res[1, ] %>%
    unlist() %>%
    unique() %>%
    round(5) %>%
    expect_identical(c(2, 13.56163, 157.94297, NA))

  res <- calculate_within_inertia(df0, cl = cl_full)
  expect_equal(dim(res), c(149, 151))

  res <- calculate_within_inertia(df, cl = cl_kmeans)
  as.numeric(res[1, 2:3]) %>%
    expect_equal(clustering_kmeans$result[[1]]$withinss)
})

# ---- TEST: DISTANCE PER VARIABLE ----
test_that("distPerVariable", {
  # Calculate distances per variable for the first partition
  res <- getDistPerVariable(df0, cl = cls[[1]])

  # Distances should be symmetric and sum to zero
  expect_equal(res[1, ], -res[2, ])
  apply(res, 2, sum) %>%
    round(7) %>%
    as.vector() %>%
    expect_equal(rep(0, 4))

  # Verify the matrix has the expected number of rows (2 clusters)
  expect_equal(nrow(res), 2)

  # Test full partition (150 clusters)
  res <- getDistPerVariable(df0, cl = cl_full[[149]])
  apply(res, 2, sum) %>%
    round(7) %>%
    as.vector() %>%
    expect_equal(rep(0, 4))
  expect_equal(nrow(res), 150)

  # Ensure consistency with scaled data
  getDistPerVariable(df, cl = cl_full[[149]]) %>%
    expect_equal(res)
})

# ---- TEST: VARIABLE CONTRIBUTIONS ----
test_that("getCtrVar", {
  # Calculate variable contributions to clustering for 3 clusters (original and scaled data)
  res <- getCtrVar2(df0, cl = cls[[2]], n_cluster = 3, scale = FALSE)
  res2 <- getCtrVar2(df, cl = cls[[2]], n_cluster = 3)

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
  res2 <- getCtrVar(df0, cl = cls[[2]], n_cluster = 3)
  expect_equal(res, res2)
  getCtrVar(df, cl = cls[[2]], n_cluster = 3) %>%
    expect_equal(res2)

  # Check specific contribution values for the first variable
  res[, 1] %>%
    round(6) %>%
    as.vector() %>%
    expect_equal(c(0.32801820, 0.02882824, 0.34358365) %>% round(6))

  # Validate results for full partition (150 clusters)
  getCtrVar2(df0, cl = cl_full[[149]], n_cluster = 150, scale = FALSE) %>%
    dim() %>%
    expect_equal(c(150, 4))
})

# ---- TEST: PARTITION-WISE DISTANCE CONTRIBUTIONS ----
test_that("getPdis", {
  # Calculate partition distance for 3 clusters
  res <- getPdis(df, cl = cls[[2]], n_cluster = 3)

  # Partition distance should match sum of variable contributions
  getCtrVar(df, cl = cls[[2]], n_cluster = 3) %>%
    apply(2, sum) %>%
    expect_equal(res)

  # Test partition distances for multiple partitions
  res2 <- getPdisPerPartition(df, cl = cls)
  expect_equal(res2[2, ], res)

  # With n = 7
  res2 <- getPdisPerPartition(df, cl = cl_full, max_cluster = 7)
  expect_equal(dim(res2), c(6, 4))

  # Specific value checks
  res2[, 1] %>%
    round(7) %>%
    as.vector() %>%
    expect_equal(c(0.4871557, 0.7004301, 0.7901667, 0.8192104, 0.8425190, 0.8660134))
})
