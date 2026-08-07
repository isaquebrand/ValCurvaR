test_that("OLS retrieves a known linear relationship", {
  d <- validar_curva(curva_homocedastica, concentracao, sinal, replica)
  fit <- ajustar_curva(d, "ols")
  expect_s3_class(fit, "valcurva_fit")
  expect_equal(unname(coef(fit$modelo)), c(0, 2), tolerance = 0.08)
  expect_equal(fit$metodo, "ols")
})

test_that("automatic fitting uses OLS for homoscedastic data", {
  d <- validar_curva(curva_homocedastica, concentracao, sinal, replica)
  fit <- ajustar_curva(d)
  expect_equal(fit$metodo, "ols")
  expect_true(is.list(fit$diagnosticos$falta_ajuste))
})

test_that("WLS uses experimental variance weights by level", {
  bruto <- data.frame(
    x = rep(1:5, each = 3),
    y = c(2, 2.1, 1.9, 4, 4.3, 3.7, 6, 6.6, 5.4, 8, 8.9, 7.1, 10, 11.2, 8.8)
  )
  d <- validar_curva(bruto, x, y)
  fit <- ajustar_curva(d, "wls", "variancia_nivel")
  expect_equal(fit$metodo, "wls")
  expect_gt(length(unique(fit$pesos)), 1)
})

test_that("validation rejects insufficient levels", {
  d <- data.frame(x = c(1, 1, 2, 2), y = c(1, 1, 2, 2))
  expect_error(validar_curva(d, x, y), "tres niveis")
})
