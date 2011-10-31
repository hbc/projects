library(testthat)

context("shRNA loading and examination")

in.data <- data.frame(shrna.id = c(rep("sh1", 3), rep("sh2", 3)),
                      d.3 = c(50, 100, 50, 5000, 5500, 3500),
                      w.3 = c(50, 75, 100, 8000, 2500, 4500),
                      accession = rep("A", 6),
                      gene.symbol = rep("", 6),
                      replicate = rep(1:3, 2))

test_that("Loading multiple shRNAs targetting a single accession", {
  reorg.data <- loadByTarget(in.data)
  expect_equal(names(reorg.data), c("accession", "spread"))
  expect_equal(as.character(reorg.data[1,1]), "A")
  expect_equal(reorg.data[1,2], 0.9866666667)
})


interactiveDevel <- function() {
  load_all("tp_shrna")
  test("tp_shrna")
}
