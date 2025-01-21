#  Estimação dos parâmetros de uma cadeia estocástica via Inferência Bayesiana

Esse projeto contém os arquivos utilizados para estimação dos parâmetros de uma cadeias estocástica via inferência Bayesiana por simulação utilizando o método MCMC Gibbs Sampling. 

São criadas duas funções:

* `GS_dirichlet`: recebe uma amostra gerada por uma cadeia estocástica, o parâmetro em interesse, o espaço de estados e o número de interações e a saída são as amostras geradas.
* `dirichlet_mean`: recebe a amostra gerada pelo `GS_dirichlet` o parâmetro em interesse e o espaço de estados e a saída é a estimação do parâmetro em interesse.

Por fim, é feito um estudo comparando as estimativas Bayesianas com as obtidas pelo método de máxima verossimilhança. 


<details>
<summary> ~\codes: </summary>

* `mcmc-dirichlet.R`: Script das funções
* `mcmc-dirichlet.Rmd`: Estudo de simualação
* `mcmc-dirichlet.html`: Relatório gerado
</details>

<details>
<summary> ~\plots: </summary>

Imagens geradas para diagnóstico das amostras geradas.

</details>
