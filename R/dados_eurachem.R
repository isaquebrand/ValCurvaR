#' Eurachem A5.2 calibration data
#'
#' Fifteen individual absorbance measurements at five cadmium calibration
#' levels. Values are transcribed from Table A5.2 of the Eurachem/CITAC
#' uncertainty guide and are supplied in long format to preserve replicates.
#'
#' @return A data frame with `concentracao_mg_L`, `replica`, and `absorbancia`.
#' @examples
#' dados <- dados_eurachem_a52()
#' curva <- validar_curva(dados, concentracao_mg_L, absorbancia, replica)
#' ajuste <- ajustar_curva(curva, metodo = "ols")
#' grafico_calibracao(ajuste)
#' @export
dados_eurachem_a52 <- function() {
  utils::read.csv(
    system.file("extdata", "eurachem_a52.csv", package = "ValCurvaR"),
    check.names = FALSE
  )
}
