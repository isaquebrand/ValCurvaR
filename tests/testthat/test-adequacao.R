test_that("adequacy assessment does not approve unspecified criteria", {
  d <- validar_curva(curva_homocedastica, concentracao, sinal, replica)
  fit <- ajustar_curva(d, "ols")
  z <- avaliar_adequacao(fit)
  expect_match(z$conclusao, "Sem decisao")
  expect_true(all(z$criterios$status == "nao avaliado"))
})

test_that("adequacy assessment applies explicit limits", {
  d <- validar_curva(curva_homocedastica, concentracao, sinal, replica)
  fit <- ajustar_curva(d, "ols")
  z <- avaliar_adequacao(fit, list(cv_max_percentual = 6, erro_relativo_retrocalculo_max_percentual = 6))
  expect_true(all(z$criterios$status[z$criterios$criterio %in% c("CV maximo por nivel (%)", "erro relativo medio maximo por nivel (%)")] == "atende"))
  audit <- relatorio_auditoria(fit, list(cv_max_percentual = 6))
  expect_true(is.matrix(audit$modelo$covariancia_coeficientes))
})
