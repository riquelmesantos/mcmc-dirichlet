# Pacotes Necessários
library(BCT)
library(LaplacesDemon)
library(ggplot2)
library(tidyr)
library(dplyr)

# ================================ FUNÇÕES AUXILIARES ============================================
### a_s(j): Número de ocorrências de uma sequêcia seguida de um simbolo específico na amostra
# s: sequencia 
# j: simbolo
# x: amostra

a_sj <- function(s, j, x) {
  n <- length(x)
  count <- 0
  for (t in (length(s) + 1):n) {
    if (identical(x[(t - length(s)):(t - 1)], s) && x[t] == j) {
      count <- count + 1
    }
  }
  return(count)
}

### $M_s$: Número de ocorrências de uma sequêcia na amostra
# s: sequencia a ser contada
# A: alfabeto
# x: amostra

M_s <- function(s, A, x) {
  count <- 0
  for (j in A) {
    count <- count + a_sj(s, j, x)
  }
  return(count)
}

### $\hat{\theta}_s(j)$: Probabilidade estimada via máxima verossimilhança
# s: sequencia
# j: simbolo
# A: alfabeto
# x: amostra

theta_sj <- function(j, s, A, amostra){
  den <- M_s(s, A, amostra)
  
  if(den == 0) return(1/length(A))
  
  num <- a_sj(s, j, amostra)
  
  return(num/den)
}

### Matriz de transição estimada via máxima verossilhança
# tau: contextos
# A: alfabeto
# x: amostra

Q_estim <- function(tau, A, x){
  Q = matrix(0, nrow = length(tau), ncol = length(A), dimnames = list(tau, A))
  for (s in 1:length(tau)) {
    for (j in 1:length(A)) {
      s_aux <- as.numeric(strsplit(tau[s], "")[[1]])
      
      den <- M_s(s_aux, A, x)
      if (den == 0){
        Q[s, j] <- 1 / length(A)
      }
      if (den > 0){
        num <- a_sj(s_aux, A[j], x)
        Q[s, j] <- num / den
      }
    }
  }
  return(Q)
}

# ================================ MCMC Dirichlet ============================================

### MCMC Dirichlet
# x: amostra
# s: sequencia
# A: alfabeto
# num_iter: número de iterações

GS_dirichlet <- function(x, s, A, num_iter = 1000) {
  A <- as.character(A)
  a_s <- sapply(A, function(j) a_sj(as.numeric(strsplit(s, "")[[1]]), j, x))
  
  t <- length(a_s)
  theta <- c()
  samples <- matrix(NA, nrow = num_iter, ncol = t)
  
  for (iter in 1:num_iter) {
    for (i in 1:(t-1)) {
      sum_a_s <- sum(a_s) - a_s[i]
      theta[i] <- rbeta(1, a_s[i], sum_a_s)
    }
    theta[t] <- 1 - sum(theta[1:t-1])
    
    samples[iter, ] <- theta
  }
  
  return(samples)
}

dirichlet_mean <- function(x, s, A){
  A <- as.character(A)
  a_s <- sapply(A, function(j) a_sj(as.numeric(strsplit(s, "")[[1]]), j, x))
  return((a_s+(1/2)) / sum(a_s + (1/2)))
}

## Exemplo 1
### Simulação da cadeia de Markov:
tau1 <- list("000" = c(0.5, 0.5), 
             "001" = c(0.4, 0.6),
             "010" = c(0.7, 0.3),
             "100" = c(0.2, 0.8),
             "011" = c(0.9, 0.1),
             "101" = c(0.3, 0.7),
             "110" = c(0.8, 0.2),
             "111" = c(0.1, 0.9))

contextos1 <- names(tau1)
A1 <- c("0", "1")
set.seed(999)
amostra1 <- as.numeric(strsplit(generate_data(tau1, 5000), "")[[1]])

orig1 <- matrix(unlist(tau1), nrow = 8, byrow = TRUE)
rownames(orig1) <- contextos1
colnames(orig1) <- A1

cat("Matriz que gerou os dados: \n")
print(orig1)

### Estimação por máxima verossimilhança dos parâmetros
estim_MLE1 <- Q_estim(contextos1, A1, amostra1)
cat("Estimação por máxima verossimilhança: \n")
print(estim_MLE1)


### Estimação via MCMC Dirichlet
num_iter <- 5000
set.seed(999)
e_000 <- GS_dirichlet(amostra1, "000", A1, num_iter)
e_001 <- GS_dirichlet(amostra1, "001", A1, num_iter)
e_010 <- GS_dirichlet(amostra1, "010", A1, num_iter)
e_100 <- GS_dirichlet(amostra1, "100", A1, num_iter)
e_011 <- GS_dirichlet(amostra1, "011", A1, num_iter)
e_101 <- GS_dirichlet(amostra1, "101", A1, num_iter)
e_110 <- GS_dirichlet(amostra1, "110", A1, num_iter)
e_111 <- GS_dirichlet(amostra1, "111", A1, num_iter)

estim_GS1 <- matrix(c(mean(e_000[,1]), mean(e_000[,2]), 
                      mean(e_001[,1]), mean(e_001[,2]),
                      mean(e_010[,1]), mean(e_010[,2]),
                      mean(e_100[,1]), mean(e_100[,2]), 
                      mean(e_011[,1]), mean(e_011[,2]), 
                      mean(e_101[,1]), mean(e_101[,2]), 
                      mean(e_110[,1]), mean(e_110[,2]), 
                      mean(e_111[,1]), mean(e_111[,2])), nrow = 8, byrow = TRUE)
rownames(estim_GS1) <- contextos1
colnames(estim_GS1) <- A1

estim_DirMean1 <- matrix(c(dirichlet_mean(amostra1, "000", A1),
                           dirichlet_mean(amostra1, "001", A1),
                           dirichlet_mean(amostra1, "010", A1),
                           dirichlet_mean(amostra1, "100", A1),
                           dirichlet_mean(amostra1, "011", A1),
                           dirichlet_mean(amostra1, "101", A1),
                           dirichlet_mean(amostra1, "110", A1),
                           dirichlet_mean(amostra1, "111", A1)), nrow = 8, byrow = TRUE)
rownames(estim_DirMean1) <- contextos1
colnames(estim_DirMean1) <- A1

### Resultados

cat("Estimação pelo GS: \n")
print(estim_GS1)
cat("Diferença entre os estimadores: \n")
print(abs(estim_GS1 - estim_MLE1))

cat("Erro quadrático médio MLE: \n")
print(unname(colMeans((estim_MLE1 - orig1)^2)[1]))

cat("Erro quadrático médio GS: \n")
print(unname(colMeans((estim_GS1 - orig1)^2)[1]))

cat("Erro quadrático médio E(X): \n")
print(unname(colMeans((estim_DirMean1 - orig1)^2)[1]))
