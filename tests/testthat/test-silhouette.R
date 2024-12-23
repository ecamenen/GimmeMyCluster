data("iris")
df0 <- select(iris, -Species)
df <- scale(df0)
dist <- dist(df)
clustering_ward <- getCAH(3, dist)
cls <- getClusterPerPart(c = clustering_ward)
cl_full <- getClusterPerPart(c = clustering_ward, n = 150)
silhouette_k <- getSilhouette(cls[[1]], dist)
silhouette_148 <- getSilhouette(cl_full[[148]], dist)
silhouettes <- getSilhouettePerPart(cl_full[-149], dist)

test_silhouette <- function(x) {
  expect_s3_class(x, "silhouette")
  min(x[, 3]) %>%
    expect_gte(-1)
  max(x[, 3]) %>%
    expect_lte(1)
  expect_length(x[, 1], 150)
}

test_that("silhouhette", {
  unique(silhouette_k[, 1]) %>%
    expect_equal(seq(2))
  head(silhouette_k[, 3]) %>%
    as.numeric() %>%
    round(7) %>%
    expect_equal(c(0.7746877, 0.7744370, 0.7739863, 0.7699909, 0.7691272, 0.7669716))
  test_silhouette(silhouette_k)
  test_silhouette(silhouette_148)
})

test_that("silhouhettePerPart", {
  expect_equal(silhouettes[[1]], silhouette_k)
  expect_equal(silhouettes[[148]], silhouette_148)
})

test_that("meanSilhouette", {
  getMeanSilhouettePerPart(silhouettes) %>%
    .[148] %>%
    expect_equal(silhouette_148[, 3] %>% mean())
})
