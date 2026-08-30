# ValCurvaR

Ferramentas transparentes para estudo de linearidade e validação de curvas
analíticas em R. O pacote preserva as replicatas, avalia falta de ajuste e
heterocedasticidade, compara OLS e WLS e produz gráficos de calibração,
resíduos, variância e incerteza de predição.

## Exemplo mínimo

```r
dados <- validar_curva(dados_brutos, concentracao, sinal, replica)
ajuste <- ajustar_curva(dados, metodo = "auto")
diagnosticos <- diagnosticar_curva(ajuste)
painel_calibracao(ajuste)
```

O resultado automático é uma recomendação baseada em diagnósticos; a decisão
de aprovar a faixa e o modelo continua documentada pelo laboratório.

## Diagnosticos adicionais

`diagnosticar_curva()` reune Shapiro-Wilk, Anderson-Darling,
Kolmogorov-Smirnov com correcao de Lilliefors e Ryan-Joiner para os residuos;
Brown-Forsythe, Breusch-Pagan, Goldfeld-Quandt e Cochran para variancia;
Durbin-Watson e Breusch-Godfrey para independencia; e residuos padronizados e
studentizados, alavancagem, Cook, DFFITS e DFBETAS para influencia.

```r
diagnosticos$normalidade
diagnosticos$independencia
grafico_qq(ajuste, semente = 2026)
grafico_influencia(ajuste)
analisar_sensibilidade(ajuste)
```

`analisar_sensibilidade()` ajusta modelos deixando uma observacao de fora
somente para medir seu impacto nos estimadores. Ela nunca altera os dados nem
recomenda a exclusao automatica de uma medicao. O pacote nao transforma os
dados para tentar obter normalidade.

## Exemplo Eurachem A5.2

```r
dados <- dados_eurachem_a52()
curva <- validar_curva(dados, concentracao_mg_L, absorbancia, replica)
ajuste_ols <- ajustar_curva(curva, metodo = "ols")
ajuste_wls <- ajustar_curva(curva, metodo = "wls")
comparar_modelos(ajuste_ols, ajuste_wls)
painel_calibracao(ajuste_wls)
```
