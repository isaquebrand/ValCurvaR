# ValCurvaR

Ferramentas transparentes para estudo de linearidade e validação de curvas
analíticas em R. O pacote preserva as replicatas, avalia falta de ajuste e
heterocedasticidade, compara OLS e WLS e produz gráficos de calibração,
resíduos, variância e incerteza de predição.

## Exemplo mínimo

```r
dados <- validar_curva(dados_brutos, concentracao, sinal, replica)
ajuste <- ajustar_curva(dados, metodo = "auto")
diagnosticar_curva(ajuste)
painel_calibracao(ajuste)
```

O resultado automático é uma recomendação baseada em diagnósticos; a decisão
de aprovar a faixa e o modelo continua documentada pelo laboratório.

## Exemplo Eurachem A5.2

```r
dados <- dados_eurachem_a52()
curva <- validar_curva(dados, concentracao_mg_L, absorbancia, replica)
ajuste_ols <- ajustar_curva(curva, metodo = "ols")
ajuste_wls <- ajustar_curva(curva, metodo = "wls")
comparar_modelos(ajuste_ols, ajuste_wls)
painel_calibracao(ajuste_wls)
```
