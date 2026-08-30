#' Generate a complete calibration-validation report in PDF
#'
#' @param ajuste A `valcurva_fit` object.
#' @param arquivo Output PDF path.
#' @param titulo Title shown in the report.
#' @param criterios Optional acceptance criteria forwarded to [avaliar_adequacao()].
#' @param quiet Logical; suppress R Markdown progress messages.
#' @return Invisibly, the normalized path to the generated PDF.
#' @export
gerar_relatorio_pdf <- function(ajuste, arquivo = "relatorio-calibracao.pdf",
                                titulo = "Relatorio de validacao da curva de calibracao",
                                criterios = list(), quiet = TRUE) {
  if (!inherits(ajuste, "valcurva_fit")) .valcurva_abort("`ajuste` deve ser um objeto `valcurva_fit`.")
  if (!requireNamespace("rmarkdown", quietly = TRUE) || !requireNamespace("knitr", quietly = TRUE))
    .valcurva_abort("Instale os pacotes `rmarkdown` e `knitr` para gerar o PDF.")
  template <- system.file("rmarkdown", "relatorio-calibracao.Rmd", package = "ValCurvaR")
  if (!nzchar(template)) .valcurva_abort("O modelo R Markdown nao foi encontrado; reinstale o pacote.")
  arquivo <- normalizePath(arquivo, winslash = "/", mustWork = FALSE)
  diretorio <- dirname(arquivo)
  if (!dir.exists(diretorio)) dir.create(diretorio, recursive = TRUE, showWarnings = FALSE)
  resultado <- rmarkdown::render(template, output_file = basename(arquivo), output_dir = diretorio,
    params = list(ajuste = ajuste, titulo = titulo, criterios = criterios),
    envir = new.env(parent = globalenv()), quiet = quiet)
  invisible(normalizePath(resultado, winslash = "/", mustWork = TRUE))
}
