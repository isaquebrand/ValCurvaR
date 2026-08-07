test_that("GUM back-calculation returns concentration and uncertainty", {
  d <- validar_curva(curva_homocedastica, concentracao, sinal, replica)
  fit <- ajustar_curva(d, "ols")
  z <- retrocalcular_concentracao(fit, c(4, 8))
  expect_equal(z$concentracao, c(2, 4), tolerance = 0.1)
  expect_true(all(z$U_expandida > 0))
  expect_true(all(is.finite(z$contribuicao_cov_intercepto_inclinacao)))
})

test_that("GUM propagation retains the fitted coefficient covariance", {
  d <- validar_curva(curva_homocedastica, concentracao, sinal, replica)
  fit <- ajustar_curva(d, "ols")
  z <- retrocalcular_concentracao(fit, 4, u_sinal = 0.1, fator_cobertura = 2)
  V <- vcov(fit$modelo); x <- z$concentracao
  expected <- sqrt((0.1 / coef(fit$modelo)[2])^2 +
    V[1, 1] / coef(fit$modelo)[2]^2 +
    (x / coef(fit$modelo)[2])^2 * V[2, 2] +
    2 * (-1 / coef(fit$modelo)[2]) * (-x / coef(fit$modelo)[2]) * V[1, 2])
  expect_equal(z$u_padrao, unname(expected))
})
