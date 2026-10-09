context("Test .unify_de_results")

## `.unify_de_results()` is the one place the four DE wrappers' shared format
## is defined, so pin the columns, the row names, and that `padj` is plain BH
## computed only over the genes with a p-value. Counting a gene the method
## did not test would inflate every other gene's adjustment.
test_that(".unify_de_results builds the shared table with BH over non-NA p-values", {
  pvalue_vec <- c(0.001, 0.01, NA, 0.5, 0.04)
  res <- .unify_de_results(gene = paste0("g", 1:5),
                           logFC = c(2, -1, NA, 0.1, 0.5),
                           se = c(0.5, 0.3, NA, 0.2, 0.2),
                           pvalue = pvalue_vec)

  expect_true(is.data.frame(res))
  expect_equal(colnames(res), c("gene", "logFC", "se", "pvalue", "padj"))
  expect_equal(rownames(res), paste0("g", 1:5))
  expect_true(is.character(res$gene))
  expect_equal(res$padj, stats::p.adjust(pvalue_vec, method = "BH"))
  expect_true(is.na(res$padj[3]))
  # four tests, not five: the NA gene does not count toward m
  expect_equal(res$padj[1], 0.001 * 4)
})

test_that(".unify_de_results rejects mismatched lengths and duplicated genes", {
  expect_error(.unify_de_results(gene = c("a", "b"),
                                 logFC = 1,
                                 se = c(1, 1),
                                 pvalue = c(0.1, 0.2)))
  expect_error(.unify_de_results(gene = c("a", "a"),
                                 logFC = c(1, 1),
                                 se = c(1, 1),
                                 pvalue = c(0.1, 0.2)))
})
