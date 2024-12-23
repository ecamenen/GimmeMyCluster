data("iris")
df0 <- select(iris, -Species)
df <- scale(df0)

# Calculate Euclidean and Manhattan distances between observations
dist_euc <- dist(df)  # Euclidean distance
dist_man <- dist(df, method = "manhattan")  # Manhattan distance

# Perform hierarchical clustering using Ward's method with Euclidean distance
clustering_ward <- getCAH(3, dist_euc)

# Extract clustering results for different partitions
cls <- getClusterPerPart(c = clustering_ward)  # Default number of partitions
cl_full <- getClusterPerPart(c = clustering_ward, n = 150)  # Full partitioning with 150 clusters

# Test clustering
test_that("clustering", {
  # Ensure the result is an hclust object
  expect_s3_class(clustering_ward, "hclust")

  # Check that Ward's method is correctly applied
  expect_identical(clustering_ward$method, "ward.D2")

  # Test Partitioning Around Medoids (PAM) clustering on distance object
  expect_s3_class(getCNH(d = dist_euc), "pam")

  # Test k-means clustering
  clustering <- getClassif(t = 2, df = df)
  expect_identical(clustering$method, "kmeans")

  # Verify that clustering results for 2 to 6 clusters are k-means objects
  sapply(2:6, function(x) class(clustering[[x]])) %>%
    unique() %>%
    expect_identical("kmeans")
})

# Function to test the structure of clustering results
test_clusters <- function(res) {
  expect_length(res, 5)  # Ensure the result contains 5 elements (for 2 to 6 clusters)
  sapply(res, max) %>%
    expect_identical(2:6)  # Verify the maximum cluster labels correspond to 2 to 6
  sapply(res, length) %>%
    unique() %>%
    expect_equal(150)  # Check that each partition contains 150 observations
}

# Test clustering results
test_that("cluster", {
  test_clusters(cls)  # Test default partitions

  # Ensure the cluster sizes match the expected counts
  table(cls[[1]]) %>%
    as.vector() %>%
    expect_equal(c(49, 101))

  # Confirm that Euclidean distance was used for the initial clustering
  expect_identical(clustering_ward$dist.method, "euclidean")

  # Test clustering with Manhattan distance
  clustering_man <- getCAH(3, dist_man)
  expect_identical(clustering_man$dist.method, "manhattan")

  # Test partitioning for Manhattan distance
  getClassif(t = 1, d = dist_man) %>%
    getClusterPerPart(c = .) %>%
    test_clusters()

  # Test partitioning for k-means
  getClassif(t = 2, df = df) %>%
    getClusterPerPart(c = .) %>%
    test_clusters()
})

# Test distance calculations
test_that("distance", {
  expect_s3_class(dist_euc, "dist")  # Ensure Euclidean distance is a dist object
  expect_no_message(checkEuclidean(dist_euc))  # Confirm it's a valid Euclidean distance

  # Verify that all non-zero distances are positive
  expect_gt(min(dist_euc[dist_euc != 0]), 0)

  # Ensure Manhattan distance is also a dist object
  expect_s3_class(dist_man, "dist")

  # Check that an error is thrown if a non-Euclidean distance is used incorrectly
  expect_error(checkEuclidean(dist_man), "Distance should be euclidean with this classification method.")
})

# Test agglomeration coefficient for clustering
test_that("coefAggl", {
  coef <- getCoefAggl(clustering_ward) %>%
    round(7)

  # Check that the coefficient is as expected and within bounds (0 < coef < 1)
  expect_identical(coef, 0.9867558)
  expect_gt(coef, 0)
  expect_lt(coef, 1)
})

# Test cluster coloring functionality
test_that("colorClusters", {
  colorClusters(cls[[2]]) %>%
    unique() %>%
    expect_identical(c("#991A80FF", "#809926FF", "#0000FFFF"))  # Ensure correct cluster colors
})

# Test cluster centroid calculation
test_that("clusterCentroids", {
  res <- getClusterCentroids(df0, cls[[1]])

  # Check that the centroid matrix has the correct dimensions (2 clusters, 4 variables)
  expect_equal(dim(res), c(2, 4))

  # Verify the centroid values for the first variable
  as.vector(res[, 1]) %>%
    round(6) %>%
    expect_identical(c(5.016327, 6.244554))

  # Test centroids for the full partition with 150 clusters
  getClusterCentroids(df0, cl_full[[149]]) %>%
    dim() %>%
    expect_equal(c(150, 4))
})
