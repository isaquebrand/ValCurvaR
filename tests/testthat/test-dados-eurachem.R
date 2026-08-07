test_that("Eurachem A5.2 data retain all individual replicates", {
  d <- dados_eurachem_a52()
  expect_equal(nrow(d), 15)
  expect_equal(length(unique(d$concentracao_mg_L)), 5)
  expect_equal(d$absorbancia[d$concentracao_mg_L == 0.9], c(0.215, 0.230, 0.216))
})
