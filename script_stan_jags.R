#Limpando os dados dos sistema
rm(list=ls(all=TRUE))

#Pacotes utilizados no trabalho
library(tidyverse)
library(rstan)
library(coda)
library(rbenchmark)
library(rjags)
library(R2jags)
library(ggmcmc)
library(MCMCvis)
library(bayesplot)
require(rjags)
require(R2jags)

#Importacao dos dados, mudar o diretorio conforme a localizacao do aquivo lasers.txt
setwd("C:/Users/guilh/OneDrive/Área de Trabalho/Pós Graduação/TCC/Projeto R TCC/TCC Jags")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

wd <- getwd()

lasers <- rio::import("lasers.txt") 
glimpse(lasers)

# no artigo aparentemente não está sendo feita essa mudança de escala,
# utilize C <- 1 para obter resultados similares ao do artigo
C <- 250 

lasers <- lasers %>%
  mutate(
    hora = hora/C
  )
glimpse(lasers)


alpha <- 0.1
Df <- 10

a_beta <- 0.01
b_beta <- 0.01
a_lambda <- 0.01
b_lambda <- 0.01
a_tau <- 0.01
b_tau <- 0.01
#Grafico para análise exploratória dos dados:

graph_laser <- ggplot(lasers, aes(x=hora, y=degradacao, group = unidade)) +
  xlab("Horas") + ylab("Degradação (%)")+
  geom_point() +
  geom_path() +
  geom_abline(intercept = Df, slope = 0, color = "blue")+
  theme_bw()+
  theme(axis.title = element_text(size = 16))
graph_laser
  
# parâmetros de interesse:
parameters <- c("beta","lambda", "sigma", "r4500", "talpha")

data <- with(
  lasers, 
  list(
    hours = hora, units = unidade, y = degradacao, n = length(degradacao), 
    L = max(unidade), alpha = alpha, Df = Df, C = C,
    a_beta = a_beta, b_beta = b_beta, 
    a_lambda = a_lambda, b_lambda = b_lambda,
    a_tau = a_tau, b_tau = b_tau
  )
)


#Carregando os arquivos do JAGS e Stan
#Uma vez carregado comente a linha 74
mod_jags <- file.path(wd, "weibull.jags")
mod_stan <- stan_model(file="weibull.stan")
#Uma vez carregados comente a linha 79
mod_stan <- readRDS("weibull.rds")
 #file.show(mod_jags)

#Execucao do software Jags

run_jags <- function(){ 
  set.seed(1234567890)
  inits <- function(){list(theta = rexp(max(lasers$unidade)), beta = rexp(1), lambda = rexp(1), tau = rexp(1))}
  output <- jags(
    data=data, 
    inits=inits, 
    parameters=parameters, 
    n.iter=2000, n.thin=1, n.burnin=1000, n.chains=4, 
    model.file=mod_jags, 
    working.directory= wd, 
    DIC=FALSE
  )  
  return(output)  
}

#Execucao do software Stan

run_stan <- function(){
  set.seed(1234567890)
  output <- sampling(mod_stan, data)
  return(output)
}

fit_jags <- run_jags()
fit_stan <- run_stan()


print(fit_jags)
print(fit_stan, parameters)

# comparando o tempo computacional replicacoes = 1:
comp1 <- benchmark(
replications = 1,
  jags = run_jags(),
  stan = run_stan()
)
# comparando o tempo computacional replicacoes = 10:
comp10 <- benchmark(
  replications = 10,
  jags = run_jags(),
  stan = run_stan()
)
# comparando o tempo computacional replicacoes = 50:
comp50 <- benchmark(
  replications = 50,
  jags = run_jags(),
  stan = run_stan()
)

# criando uma função generica para extrair o tamamnho efetivo da cadeia:
effective_size <- function(object, parameters, ...) UseMethod("effective_size")

effective_size.rjags <- function(object, parameters){
  n_eff <- object$BUGSoutput$summary[parameters, "n.eff"]
  return(n_eff)
}

effective_size.stanfit <- function(object, parameters){
  tab <- rstan::summary(object)
  n_eff <- tab$summary[parameters, "n_eff"]
  return(n_eff)
}


neff_jags<- effective_size(fit_jags, parameters)
neff_stan <- effective_size(fit_stan, parameters)

# comparando tamanho efetivo da cadeia:
neff_stan/neff_jags


#Calculo do HPD

H1 <- HPDinterval(samp_jags)
H2 <- HPDinterval(samp_stan)
Ap1 <- apply(H1,  1, diff)
Ap2 <- apply(H2,  1, diff)
Ap1/Ap2

##########
#Amostras a posteriori 
samp_stan <- as.mcmc(with(extract(fit_stan), cbind(beta, lambda, sigma, r4500, talpha)))
samp_jags <- as.mcmc(as.data.frame(fit_jags$BUGSoutput$sims.list))[, parameters]


#GRAFICOS
#funcao extract para extrair dados do Stanfit
stan_data <- extract(fit_stan,permuted = TRUE)

#funcao MCMCchains para extrair os dados do Jags
jags_data <- as.data.frame(MCMCchains(fit_jags))

#data frame dos graficos STAN e JAGS
graph_stan <- data.frame(
  k = 1:length(stan_data$beta),
  beta = stan_data$beta,
  lambda = stan_data$lambda,
  sigma = stan_data$sigma,
  software = "stan"
)

graph_jags <- data.frame(
  k = 1:length(jags_data$beta),
  beta = jags_data$beta,
  lambda = jags_data$lambda,
  sigma = jags_data$sigma,
  software = "jags"
)

#Criando um data set com todos os dados

amostra_graph <- bind_rows(graph_jags,graph_stan)
glimpse(amostra_graph)


amostra_graph_long <- amostra_graph %>%
  pivot_longer(
    cols = c("beta", "lambda", "sigma"),
    names_to = "parametro",
    values_to = "amostra_graph"
  )
glimpse(amostra_graph_long)

ggplot(amostra_graph_long, aes(x = k, y = amostra_graph, color = software)) +
  geom_path(alpha = 0.7) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  facet_grid(parametro~.,scales = "free_y",
             switch = "y",
             labeller = as_labeller(c(beta = "Beta",lambda = "Lambda", sigma = "sigma")))+
  ylab(NULL) +
  theme(strip.background = element_blank(), # remove the background
        strip.placement = "outside") # put labels to the left of the axis text

#Obtendo dados das cadeias de markov de cada modelo

stan_chains <- as.data.frame(extract(fit_stan,permuted = FALSE))

graph_stan_chains_beta <- data.frame(
  j = 1:length(stan_chains$`chain:1.beta`),
  uper_parameter = rep("beta",1000),
  chain1 = stan_chains$`chain:1.beta`,
  chain2 = stan_chains$`chain:2.beta`,
  chain3 = stan_chains$`chain:3.beta`,
  chain4 = stan_chains$`chain:4.beta`,
  software = "stan"
  )
glimpse(graph_stan_chains_beta)

graph_stan_chains_lambda <- data.frame(
  j = 1:length(stan_chains$`chain:1.lambda`),
  uper_parameter = rep("lambda",1000),
  chain1 = stan_chains$`chain:1.lambda`,
  chain2 = stan_chains$`chain:2.lambda`,
  chain3 = stan_chains$`chain:3.lambda`,
  chain4 = stan_chains$`chain:4.lambda`,
  software = "stan"
)

glimpse(graph_stan_chains_lambda)
  
graph_stan_chains_sigma <- data.frame(  
  j = 1:length(stan_chains$`chain:1.lambda`),
  uper_parameter = rep("sigma",1000),
  chain1 = stan_chains$`chain:1.sigma`,
  chain2 = stan_chains$`chain:2.sigma`,
  chain3 = stan_chains$`chain:3.sigma`,
  chain4 = stan_chains$`chain:4.sigma`,
  software = "stan"
  )

glimpse(graph_stan_chains_sigma)


graph_jags_chains_beta <- data.frame(
  j = 1:length(MCMCchains(fit_jags, params = c("beta"),mcmc.list = FALSE, chain_num = 1)),
  uper_parameter = rep("beta",1000),
  chain1 = as.vector(MCMCchains(fit_jags, params = c("beta"),mcmc.list = FALSE, chain_num = 1)),
  chain2 = as.vector(MCMCchains(fit_jags, params = c("beta"),mcmc.list = FALSE, chain_num = 2)),
  chain3 = as.vector(MCMCchains(fit_jags, params = c("beta"),mcmc.list = FALSE, chain_num = 3)),
  chain4 = as.vector(MCMCchains(fit_jags, params = c("beta"),mcmc.list = FALSE, chain_num = 4)),
  software = "jags"
)
glimpse(graph_jags_chains_beta)

  
graph_jags_chains_lambda <- data.frame(  
  j = 1:length(MCMCchains(fit_jags, params = c("beta"),mcmc.list = FALSE, chain_num = 1)),
  uper_parameter = rep("lambda",1000),
  chain1 = as.vector(MCMCchains(fit_jags, params = c("lambda"),mcmc.list = FALSE, chain_num = 1)),
  chain2 = as.vector(MCMCchains(fit_jags, params = c("lambda"),mcmc.list = FALSE, chain_num = 2)),
  chain3 = as.vector(MCMCchains(fit_jags, params = c("lambda"),mcmc.list = FALSE, chain_num = 3)),
  chain4 = as.vector(MCMCchains(fit_jags, params = c("lambda"),mcmc.list = FALSE, chain_num = 4)),
  software = "jags"
)
glimpse(graph_jags_chains_lambda)  

  graph_jags_chains_sigma <- data.frame(  
  j = 1:length(MCMCchains(fit_jags, params = c("beta"),mcmc.list = FALSE, chain_num = 1)),
  uper_parameter = rep("sigma",1000),
  chain1 = as.vector(MCMCchains(fit_jags, params = c("sigma"),mcmc.list = FALSE, chain_num = 1)),
  chain2 = as.vector(MCMCchains(fit_jags, params = c("sigma"),mcmc.list = FALSE, chain_num = 2)),
  chain3 = as.vector(MCMCchains(fit_jags, params = c("sigma"),mcmc.list = FALSE, chain_num = 3)),
  chain4 = as.vector(MCMCchains(fit_jags, params = c("sigma"),mcmc.list = FALSE, chain_num = 4)),
  software = "jags"
  )

glimpse(graph_jags_chains_sigma)
  
chains_graph <- bind_rows(graph_stan_chains_beta,graph_stan_chains_lambda, graph_stan_chains_sigma,graph_jags_chains_beta,graph_jags_chains_lambda,graph_jags_chains_sigma)
glimpse(chains_graph)

colnames(chains_graph) <- c('Iteração','upper_parameter', 'cadeia 1', 'cadeia 2', 'cadeia 3', 'cadeia 4', 'software')

chains_graph_long <- chains_graph %>%
  pivot_longer(
    cols = c("cadeia 1","cadeia 2","cadeia 3","cadeia 4"),
    names_to = "Cadeias",
    values_to = "chains_graph"
  )
glimpse(chains_graph_long)


#GRAFICO DAS CADEIAS
ggplot(chains_graph_long, aes(x = Iteração, y = chains_graph, color = Cadeias)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_path(alpha = 0.5) +
  facet_grid(upper_parameter~software~.,scales = "free_y", 
             switch = "y",
             labeller = as_labeller(c(beta = "Beta",lambda = "Lambda", sigma = "Sigma", stan = "Stan", jags = "JAGS")))+
  ylab(NULL) + 
  theme(strip.background = element_blank(), # remove the background
        strip.placement = "outside")


#Grafico de Densidade:

ggplot(chains_graph_long, aes(x=chains_graph, col = Cadeias)) +  
  geom_density()+
  facet_wrap(software~upper_parameter~., scales = "free")+
  labs(y = "Densidade", x = "Valor dos Parâmetros")

#############################################
#ANALISE DE CONFIABILIDADE
#Stan
graph_stan_chains_r4500 <- data.frame(
  j = 1:length(stan_chains$`chain:1.r4500`),
  uper_parameter = rep("r4500",1000),
  chain1 = stan_chains$`chain:1.r4500`,
  chain2 = stan_chains$`chain:2.r4500`,
  chain3 = stan_chains$`chain:3.r4500`,
  chain4 = stan_chains$`chain:4.r4500`,
  software = "stan"
)
glimpse(graph_stan_chains_r4500)

graph_stan_chains_talpha <- data.frame(
  j = 1:length(stan_chains$`chain:1.talpha`),
  uper_parameter = rep("talpha",1000),
  chain1 = stan_chains$`chain:1.talpha`,
  chain2 = stan_chains$`chain:2.talpha`,
  chain3 = stan_chains$`chain:3.talpha`,
  chain4 = stan_chains$`chain:4.talpha`,
  software = "stan"
)
glimpse(graph_stan_chains_talpha)

#JAGS
graph_jags_chains_r4500 <- data.frame(
  j = 1:length(MCMCchains(fit_jags, params = c("r4500"),mcmc.list = FALSE, chain_num = 1)),
  uper_parameter = rep("r4500",1000),
  chain1 = as.vector(MCMCchains(fit_jags, params = c("r4500"),mcmc.list = FALSE, chain_num = 1)),
  chain2 = as.vector(MCMCchains(fit_jags, params = c("r4500"),mcmc.list = FALSE, chain_num = 2)),
  chain3 = as.vector(MCMCchains(fit_jags, params = c("r4500"),mcmc.list = FALSE, chain_num = 3)),
  chain4 = as.vector(MCMCchains(fit_jags, params = c("r4500"),mcmc.list = FALSE, chain_num = 4)),
  software = "jags"
)
glimpse(graph_jags_chains_r4500)

graph_jags_chains_talpha <- data.frame(
  j = 1:length(MCMCchains(fit_jags, params = c("talpha"),mcmc.list = FALSE, chain_num = 1)),
  uper_parameter = rep("talpha",1000),
  chain1 = as.vector(MCMCchains(fit_jags, params = c("talpha"),mcmc.list = FALSE, chain_num = 1)),
  chain2 = as.vector(MCMCchains(fit_jags, params = c("talpha"),mcmc.list = FALSE, chain_num = 2)),
  chain3 = as.vector(MCMCchains(fit_jags, params = c("talpha"),mcmc.list = FALSE, chain_num = 3)),
  chain4 = as.vector(MCMCchains(fit_jags, params = c("talpha"),mcmc.list = FALSE, chain_num = 4)),
  software = "jags"
)
glimpse(graph_jags_chains_talpha)

#Agrupando os dados

relia_graph <- bind_rows(graph_stan_chains_r4500,graph_stan_chains_talpha,graph_jags_chains_r4500,graph_jags_chains_talpha)
colnames(relia_graph) <- c('Iteração','upper_parameter', 'cadeia 1', 'cadeia 2', 'cadeia 3', 'cadeia 4', 'software')

#Versao long_data
relia_graph_long <- relia_graph %>%
  pivot_longer(
    cols = c("cadeia 1","cadeia 2","cadeia 3","cadeia 4"),
    names_to = "Cadeias",
    values_to = "relia_graph"
  )
glimpse(relia_graph_long)

#Traceplot

ggplot(relia_graph_long, aes(x = Iteração, y = relia_graph, color = Cadeias)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_path(alpha = 0.5) +
  facet_grid(upper_parameter~software~.,scales = "free_y", 
             switch = "y",
             labeller = as_labeller(c(r4500 = "R(4500)",talpha = "t_alpha", stan = "Stan", jags = "JAGS")))+
  ylab(NULL) + 
  theme(strip.background = element_blank(), # remove the background
        strip.placement = "outside")

#Grafico de densidade

ggplot(relia_graph_long, aes(x=relia_graph, col = Cadeias)) +  
  geom_density()+
  facet_wrap(software~upper_parameter, scales = "free",
             labeller = as_labeller(c(r4500 = "R(4500)",talpha = "t_alpha", stan = "Stan", jags = "JAGS")))+
  labs(y = "Densidade", x = "Valor dos Parâmetros")


#Analise de autocorrelacao

#JAGS
#Beta
acf(graph_jags_chains_beta$chain1, main = "JAGS - Beta Cadeia 1", cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
acf(graph_jags_chains_beta$chain2, main = "JAGS - Beta Cadeia 2")
acf(graph_jags_chains_beta$chain3, main = "JAGS - Beta Cadeia 3")
acf(graph_jags_chains_beta$chain4, main = "JAGS - Beta Cadeia 4")

#Lambda
acf(graph_jags_chains_lambda$chain1, main = "JAGS - Lambda Cadeia 1", cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
acf(graph_jags_chains_lambda$chain2, main = "JAGS - Lambda Cadeia 2")
acf(graph_jags_chains_lambda$chain3, main = "JAGS - Lambda Cadeia 3")
acf(graph_jags_chains_lambda$chain4, main = "JAGS - Lambda Cadeia 4")

#Sigma
acf(graph_jags_chains_sigma$chain1, main = "JAGS - Sigma Cadeia 1", cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
acf(graph_jags_chains_sigma$chain2, main = "JAGS - Sigma Cadeia 2")
acf(graph_jags_chains_sigma$chain3, main = "JAGS - Sigma Cadeia 3")
acf(graph_jags_chains_sigma$chain4, main = "JAGS - Sigma Cadeia 4")

#R(4500)
acf(graph_jags_chains_r4500$chain1, main = "JAGS - R(4500) Cadeia 1", cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
acf(graph_jags_chains_r4500$chain2, main = "JAGS - R(4500) Cadeia 2")
acf(graph_jags_chains_r4500$chain3, main = "JAGS - R(4500) Cadeia 3")
acf(graph_jags_chains_r4500$chain4, main = "JAGS - R(4500) Cadeia 4")

#Talpha
acf(graph_jags_chains_talpha$chain1, main = "JAGS - t_alpha Cadeia 1", cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
acf(graph_jags_chains_talpha$chain2, main = "JAGS - t_alpha Cadeia 2")
acf(graph_jags_chains_talpha$chain3, main = "JAGS - t_alpha Cadeia 3")
acf(graph_jags_chains_talpha$chain4, main = "JAGS - t_alpha Cadeia 4")

#################################################################
#Stan
#Beta
acf(graph_stan_chains_beta$chain1, main = "Stan - Beta Cadeia 1", cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
acf(graph_stan_chains_beta$chain2, main = "Stan - Beta Cadeia 2")
acf(graph_stan_chains_beta$chain3, main = "Stan - Beta Cadeia 3")
acf(graph_stan_chains_beta$chain4, main = "Stan - Beta Cadeia 4")

#Lambda
acf(graph_stan_chains_lambda$chain1, main = "Stan - Lambda Cadeia 1", cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
acf(graph_stan_chains_lambda$chain2, main = "Stan - Lambda Cadeia 2")
acf(graph_stan_chains_lambda$chain3, main = "Stan - Lambda Cadeia 3")
acf(graph_stan_chains_lambda$chain4, main = "Stan - Lambda Cadeia 4")

#Sigma
acf(graph_stan_chains_sigma$chain1, main = "Stan - Sigma Cadeia 1", cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
acf(graph_stan_chains_sigma$chain2, main = "Stan - Sigma Cadeia 2")
acf(graph_stan_chains_sigma$chain3, main = "Stan - Sigma Cadeia 3")
acf(graph_stan_chains_sigma$chain4, main = "Stan - Sigma Cadeia 4")

#R(4500)
acf(graph_stan_chains_r4500$chain1, main = "Stan - R(4500) Cadeia 1", cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
acf(graph_stan_chains_r4500$chain2, main = "Stan - R(4500) Cadeia 2")
acf(graph_stan_chains_r4500$chain3, main = "Stan - R(4500) Cadeia 3")
acf(graph_stan_chains_r4500$chain4, main = "Stan - R(4500) Cadeia 4")

#Talpha
acf(graph_stan_chains_talpha$chain1, main = "Stan - t_alpha Cadeia 1", cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
acf(graph_stan_chains_talpha$chain2, main = "Stan - t_alpha Cadeia 2")
acf(graph_stan_chains_talpha$chain3, main = "Stan - t_alpha Cadeia 3")
acf(graph_stan_chains_talpha$chain4, main = "Stan - t_alpha Cadeia 4")

####################################################################################
#Grafico de confiabilidade:

rel_stan_beta <- c(graph_stan_chains_beta$chain1,graph_stan_chains_beta$chain2,graph_stan_chains_beta$chain3,graph_stan_chains_beta$chain4)
glimpse(rel_stan_beta)

rel_stan_lambda <- c(graph_stan_chains_lambda$chain1,graph_stan_chains_lambda$chain2,graph_stan_chains_lambda$chain3,graph_stan_chains_lambda$chain4)
glimpse(rel_stan_lambda)

rel_jags_beta <- c(graph_jags_chains_beta$chain1,graph_jags_chains_beta$chain2,graph_jags_chains_beta$chain3,graph_jags_chains_beta$chain4)
glimpse(rel_jags_beta)

rel_jags_lambda <- c(graph_jags_chains_lambda$chain1,graph_jags_chains_lambda$chain2,graph_jags_chains_lambda$chain3,graph_jags_chains_lambda$chain4)
glimpse(rel_jags_lambda)


test <- extract(fit_stan, pars = c("beta", "lambda")) %>%
  as.data.frame()

t <- seq(1,50)
t

confiabilidadestan <- array(dim=c(4000,50))
confiabilidadejags <- array(dim=c(4000,50))

data_conf_jags <- rep(0,50)
data_conf_jags.025 <- rep(0,50)
data_conf_jags.975 <- rep(0,50)

data_conf_stan <- rep(0,50)
data_conf_stan.025 <- rep(0,50)
data_conf_stan.975 <- rep(0,50)

 

for (x in 1:50) {
  for (y in 1:4000) {
        confiabilidadejags[y,x]  <- (exp(-((rel_jags_lambda[y])/(Df^rel_jags_beta[y]))*((x)^rel_jags_beta[y])))
        confiabilidadestan[y,x]  <- (exp(-((rel_stan_lambda[y])/(Df^rel_stan_beta[y]))*((x)^rel_stan_beta[y])))
  }
}
glimpse(confiabilidadejags)

glimpse(confiabilidadestan)

for (i in 1:50) {
  data_conf_jags[i] <- median(confiabilidadejags[,i])
  data_conf_jags.025[i] <- quantile(confiabilidadejags[,i],0.025)
  data_conf_jags.975[i] <- quantile(confiabilidadejags[,i],0.975)
}

glimpse(data_conf_jags)
glimpse(data_conf_jags.025)
glimpse(data_conf_jags.975)

for (i in 1:50) {
  data_conf_stan[i] <- median(confiabilidadestan[,i])
  data_conf_stan.025[i] <- quantile(confiabilidadestan[,i],0.025)
  data_conf_stan.975[i] <- quantile(confiabilidadestan[,i],0.975)
}
glimpse(data_conf_stan)
glimpse(data_conf_stan.025)
glimpse(data_conf_stan.975)

tempo_correto <- t*250
tempo_correto

graf_conf_jags <- data.frame(tempo_correto,data_conf_jags,data_conf_jags.025,data_conf_jags.975)
glimpse(graf_conf_jags)

graf_conf_stan <- data.frame(tempo_correto,data_conf_stan,data_conf_stan.025,data_conf_stan.975)
glimpse(graf_conf_stan)


ggplot()+
  geom_line(data = graf_conf_jags, aes(x=tempo_correto,y=data_conf_jags), color = 'blue')+
  geom_line(data = graf_conf_jags, aes(x=tempo_correto,y=data_conf_jags.025), color = 'black', linetype="dotted")+
  geom_line(data = graf_conf_jags, aes(x=tempo_correto,y=data_conf_jags.975), color = 'black', linetype="dotted")+
  labs(y = "Confiabilidade JAGS", x = "Horas")+
  geom_abline(intercept = 0.9, slope = 0, color = "red")+
  geom_vline(intercept = 4500, xintercept = 4500, color = "green")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  scale_x_continuous( limits=c(3000, 7000),breaks=seq(3000,7000,500))

ggplot()+
  geom_line(data = graf_conf_stan, aes(x=tempo_correto,y=data_conf_stan), color = 'blue')+
  geom_line(data = graf_conf_stan, aes(x=tempo_correto,y=data_conf_stan.025), color = 'black', linetype="dotted")+
  geom_line(data = graf_conf_stan, aes(x=tempo_correto,y=data_conf_stan.975), color = 'black', linetype="dotted")+
  labs(y = "Confiabilidade Stan", x = "Horas")+
  geom_abline(intercept = 0.9, slope = 0, color = "red")+
  geom_vline(intercept = 4500, xintercept = 4500, color = "green")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  scale_x_continuous( limits=c(3000, 7000),breaks=seq(3000,7000,500))
rm(R_t)


