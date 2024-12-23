data("iris")
df0 <- select(iris, -Species)
df <- scale(df0)
dist <- dist(df)
dist2 <- dist(df, method = "manhattan")
clustering_ward <- getCAH(3, dist)
cls <- getClusterPerPart(c = clustering_ward)
cl_full <- getClusterPerPart(c = clustering_ward, n = 150)

test_that("clustering", {
  expect_s3_class(clustering_ward, "hclust")
  expect_identical(clustering_ward$method, "ward.D2")
  expect_s3_class(getCNH(d = dist), "pam")
  clustering <- getClassif(t = 2, df = df)
  expect_identical(clustering$method, "kmeans")
  sapply(2:6, function(x) class(clustering[[x]])) %>%
    unique() %>%
    expect_identical("kmeans")
})

test_clusters <- function(res) {
  expect_length(res, 5)
  sapply(res, max) %>%
    expect_identical(2:6)
  sapply(res, length) %>%
    unique() %>%
    expect_equal(150)
}

test_that("cluster cah", {
  test_clusters(cls)
  table(cls[[1]]) %>%
    as.vector() %>%
    expect_equal(c(49, 101))
  expect_identical(clustering_ward$dist.method, "euclidean")
  clustering <- getCAH(3, dist2)
  expect_identical(clustering$dist.method, "manhattan")
  getClassif(t = 1, d = dist2) %>%
    getClusterPerPart(c = .) %>%
    test_clusters()
  getClassif(t = 2, df = df) %>%
    getClusterPerPart(c = .) %>%
    test_clusters()
})

test_that("distance", {
  expect_s3_class(dist, "dist")
  checkEuclidean(dist) %>%
    expect_no_message()
  expect_gt(min(dist[dist != 0]), 0)
  expect_s3_class(dist2, "dist")
  checkEuclidean(dist2) %>%
    expect_error("Distance should be euclidean with this classification method.")
})

test_that("coefAggl", {
  res <- getCoefAggl(clustering_ward) %>%
    round(7)
  expect_identical(res, 0.9867558)
  expect_gt(res, 0)
  expect_lt(res, 1)
})

test_that("colorClusters", {
  colorClusters(cls[[2]]) %>%
    unique() %>%
    expect_identical(c("#991A80FF", "#809926FF", "#0000FFFF"))
})

test_that("clusterCentroids", {
  res <- getClusterCentroids(df0, cls[[1]])
  expect_equal(dim(res), c(2, 4))
  as.vector(res[, 1]) %>%
    round(6) %>%
    expect_identical(c(5.016327, 6.244554))
  getClusterCentroids(df0, cl_full[[149]]) %>%
    dim() %>%
    expect_equal(c(150, 4))
})
