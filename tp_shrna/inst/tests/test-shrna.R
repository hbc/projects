library(testthat)

context("shRNA loading and examination")

in_data <- data.frame(shrna.id = c(rep("sh1", 3), rep("sh2", 3)),
                      d.3 = c(50, 100, 50, 5000, 5500, 3500),
                      w.3 = c(50, 75, 100, 8000, 2500, 4500),
                      accession = rep("A", 6),
                      gene.symbol = rep("", 6),
                      replicate = rep(1:3, 2))
config <- list(min_count=0,
               model=data.frame(conditions = c(rep("d.3", 3), rep("w.3", 3)),
                                type = rep(c("standard", "suspect", "standard"), 2)))

test_that("Prepare input data for DESeq, by gene", {
  prep_data <- prepareByAccession(in_data, config)
  expect_equal(length(prep_data$model$conditions), 6)
  expect_equal(rownames(prep_data$counts), c("A"))
  expect_equal(nrow(prep_data$counts), 1)
  expect_equal(prep_data$counts[1,1], 5050)
})

test_that("Prepare input data for DESeq, by target", {
  data <- prepareByTarget(in_data, config)
  expect_equal(ncol(data$counts), 6)
  expect_equal(rownames(data$counts), c("sh1", "sh2"))
  expect_equal(data$counts[1,2], 100)
})

interactiveDevel <- function() {
  load_all("tp_shrna")
  test("tp_shrna")
}

test_that("Filter input data by counts", {
  f_data <- filterDfByCounts(in_data, 500)
  expect_equal(nrow(f_data), 3)
  expect_equal(as.character(f_data$shrna.id[1]),  "sh2")
  f_data <- filterDfByCounts(in_data, 10)
  expect_equal(nrow(f_data), 6)
})

test_that("Loading multiple shRNAs targetting a single accession", {
  reorg.data <- loadByTarget(in_data)
  expect_equal(names(reorg.data), c("accession", "spread"))
  expect_equal(as.character(reorg.data[1,1]), "A")
  expect_equal(reorg.data[1,2], 0.9866666667)
})
