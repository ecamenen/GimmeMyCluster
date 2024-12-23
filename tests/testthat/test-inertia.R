data("iris")
df0 <- select(iris, -Species)
df <- scale(df0)
clustering_ward <- dist(df) %>% getCAH(3, .)
cls <- getClusterPerPart(c = clustering_ward)
cl_full <- getClusterPerPart(c = clustering_ward, n = 150)

test_that("betweenPerPart", {
  res <- getRelativeBetweenPerPart(d = df, cl = cls)
  round(res, 5) %>%
    expect_identical(c(61.87981, 75.18729, 80.52696, 84.16555, 85.66615))
  # getRelativeBetweenPerPart(d = df0, cl = cls) %>%
  #   expect_equal(res)
  getBetweenDifferences(res) %>%
    round(6) %>%
    expect_identical(c(61.879806, 13.307485,  5.339665,  3.638599,  1.500597))
  cl_full %>%
    getRelativeBetweenPerPart(d = df, cl = ., n = 150) %>%
    .[149] %>%
    expect_identical(100)
})

test_that("betweenPerPart", {
  res <- getRelativeWithinPerCluster(cls = cls, d = df0)
  getRelativeWithinPerCluster(cls = cls, d = df) %>%
    expect_equal(res)
  apply(res, 1, function(x) sum(x, na.rm = TRUE)) %>%
    round(6) %>%
    unique() %>%
    expect_equal(1)
  expect_equal(dim(res), c(5, 6))
  res[1, ] %>%
    unique() %>%
    round(7) %>%
    expect_identical(c(0.6733333, 0.3266667, NA))
  # TODO: after getRelativeWithinPerCluster optimization, test for n = 149
})

test_that("distPerVariable", {
  res <- getDistPerVariable(d = df0, cl = cls[[1]])
  expect_equal(res[1, ], -res[2, ])
  apply(res, 2, sum) %>%
    round(7) %>%
    as.vector() %>%
    expect_equal(rep(0, 4))
  expect_equal(nrow(res), 2)
  res <- getDistPerVariable(d = df0, cl = cl_full[[149]])
  apply(res, 2, sum) %>%
    round(7) %>%
    as.vector() %>%
    expect_equal(rep(0, 4))
  expect_equal(nrow(res), 150)
  getDistPerVariable(d = df, cl = cl_full[[149]]) %>%
    expect_equal(res)
})

test_that("getCtrVar", {
  res <- getCtrVar2(cl = cls[[2]], d = df0, k = 3, scale = FALSE)
  res2 <- getCtrVar2(cl = cls[[2]], d = df, k = 3)
  unlist(res2) %>%
    as.numeric() %>%
    round(7) %>%
    setdiff(unlist(res) %>% as.numeric() %>% round(7)) %>%
    length() %>%
    expect_gt(0)
  apply(res2, 2, sum) %>%
    as.numeric() %>%
    round(6) %>%
    expect_equal(rep(1, 4))
  res2 <- getCtrVar(cl = cls[[2]], d = df0, k = 3)
  expect_equal(res, res2)
  getCtrVar(cl = cls[[2]], d = df, k = 3) %>%
    expect_equal(res2)
  res[, 1] %>%
    as.numeric() %>%
    round(6) %>%
    expect_equal(c(0.32801820, 0.02882824, 0.34358365) %>% round(6))
  getCtrVar2(cl = cl_full[[149]], d = df0, k = 150, scale = FALSE) %>%
    dim() %>%
    expect_equal(c(150, 4))
})

test_that("getCtrVar", {
  res <- getPdis(cl = cls[[2]], d = df, k = 3)
  getCtrVar(cl = cls[[2]], d = df, k = 3) %>%
    apply(2, sum) %>%
    expect_equal(res)
  res2 <- getPdisPerPartition(cls = cls, d = df)
  expect_equal(res2[2, ], res)
  res2 <- getPdisPerPartition(cls = cl_full, d = df, n = 7)
  res2[, 1] %>%
    as.numeric() %>%
    round(7) %>%
    expect_equal(c(0.4871557, 0.7004301, 0.7901667, 0.8192104, 0.8425190, 0.8660134))
  expect_equal(dim(res2), c(6, 4))
})
