# ---- DATA PREPARATION ----
data("iris")
df0 <- select(iris, -Species)
df <- scale(df0)

# Calculate Euclidean and Manhattan distances between observations
dist_euc <- dist(df)  # Euclidean distance
dist_man <- dist(df, method = "manhattan")  # Manhattan distance

# Perform hierarchical clustering using Ward's method with Euclidean distance
clustering_ward <- run_ahc(dist_euc)

# Extract clustering results for different partitions
cls <- getClusterPerPart(clustering_ward)  # Default number of partitions
cl_full <- getClusterPerPart(clustering_ward, max_cluster = 150)  # Full partition with 150 clusters

# ---- CLUSTERING TESTS ----
test_that("clustering", {
  # Verify the clustering object is of class 'hclust'
  expect_s3_class(clustering_ward, "hclust")

  # Confirm the clustering method is 'ward.D2'
  expect_identical(clustering_ward$method, "ward.D2")

  # Test Partitioning Around Medoids (PAM) with Euclidean distance
  expect_s3_class(run_nhc(dist = dist_euc), "pam")

  # Test k-means clustering
  clustering <- run_clustering(data = df, method = 2)
  expect_identical(clustering$method, "kmeans")

  # Ensure clustering results for 2 to 6 clusters are all k-means objects
  sapply(seq(5), function(x) class(clustering$results[[x]])) %>%
    unique() %>%
    expect_identical("kmeans")
})

# ---- CLUSTER STRUCTURE TESTING FUNCTION ----
test_clusters <- function(res) {
  expect_length(res, 5)  # Ensure the result contains 5 elements (for 2 to 6 clusters)
  sapply(res, max) %>% expect_identical(2:6)  # Max cluster labels should match 2 to 6
  sapply(res, length) %>% unique() %>% expect_equal(150)  # Each partition should contain 150 obs
}

# ---- TEST PARTITION RESULTS ----
test_that("cluster", {
  # Test default partition structure
  test_clusters(cls)

  # Verify specific cluster sizes (for 2 clusters)
  table(cls[[1]]) %>% as.vector() %>% expect_equal(c(49, 101))

  # Confirm that Euclidean distance was used for the initial clustering
  expect_identical(clustering_ward$dist.method, "euclidean")

  # Test clustering with Manhattan distance
  clustering_man <- run_ahc(dist_man)
  expect_identical(clustering_man$dist.method, "manhattan")

  # Test PAM clustering with Manhattan distance
  run_clustering(dist = dist_man, method = 1) %>%
    extract_clusters() %>%
    pluck("clusters") %>%
    test_clusters()

  # Test k-means clustering with standardized data
  run_clustering(data = df, method = 2) %>%
    extract_clusters() %>%
    pluck("clusters") %>%
    test_clusters()
})

# ---- DISTANCE TESTS ----
test_that("distance", {
  # Verify Euclidean distance is of class 'dist'
  expect_s3_class(dist_euc, "dist")

  # Confirm no issues with Euclidean distance
  expect_no_message(checkEuclidean(dist_euc))

  # Ensure no zero distances (except diagonal) and all values are positive
  expect_gt(min(dist_euc[dist_euc != 0]), 0)

  # Verify Manhattan distance is of class 'dist'
  expect_s3_class(dist_man, "dist")

  # Check error message for using non-Euclidean distance where Euclidean is required
  expect_error(checkEuclidean(dist_man), "Distance should be euclidean with this classification method.")
})

# ---- AGGLOMERATION COEFFICIENT TEST ----
test_that("coefAggl", {
  coef <- round(getCoefAggl(clustering_ward), 7)

  # Ensure the coefficient is within the range (0, 1)
  expect_identical(coef, 0.9867558)
  expect_gt(coef, 0)
  expect_lt(coef, 1)
})

# ---- CLUSTER COLORING TEST ----
test_that("color_cluster", {
  # Validate cluster coloring is applied correctly
  color_cluster(cls[[2]]) %>%
    unique() %>%
    expect_identical(brewer.pal(9, "Set1")[seq(3)])
})

# ---- CLUSTER CENTROID TEST ----
test_that("clusterCentroids", {
  # Compute cluster centroids for 2 clusters
  res <- calculate_centroids(df0, cls[[1]])

  # Confirm centroid matrix dimensions (2 clusters, 4 variables + rownames)
  expect_equal(dim(res), c(2, 5))
  as.vector(pull(res, 1)) %>%
    expect_identical(seq(2))

  # Check centroid values for the first variable
  as.vector(pull(res, 2)) %>%
    round(6) %>%
    expect_identical(c(5.016327, 6.244554))

  # Test centroids for 150-cluster partition
  calculate_centroids(df0, cl_full[[149]]) %>%
    dim() %>%
    expect_equal(c(150, 5))
})
