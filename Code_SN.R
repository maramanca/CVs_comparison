library("coda")

########################## DATA ################################################################
### Data is available from the public genomics data repository Gene Expression Omnibus (GEO) ###
### and accessible through GEO Series accession number GSE18392.                             ### 
################################################################################################

#### Full conditional update for w ####
update_w <- function(beta, rec_eta2, b, d) {
  
  sp <- (d+1)/2
  rt <- (1/2)*(((beta^2*rec_eta2)/b)+d)
  
  rgamma(n = 1, shape = sp, rate = rt)
}

#### Full conditional update for u ####
update_u <- function(y,  csi, beta, rec_eta2, n_oss) {
  
  m <- ((y-csi)*beta)/(1 + rec_eta2*beta^2)
  sd <- sqrt((1 + rec_eta2*beta^2)^(-1))
  
  rnorm(n = n_oss, mean = m, sd = sd)
}

#### Full conditional update for csi ####
update_csi <- function(beta, rec_eta2, u, y, n_oss) {
  
  m <- sum(y-beta*u)/n_oss
  sd <- sqrt((n_oss*rec_eta2)^(-1))
  
  rnorm(n = 1, mean = m, sd = sd)
}

#### Full conditional update for beta ####
update_beta <- function(csi, rec_eta2, w, u, y, b) {
  
  m <- sum((y-csi)*u)/((w/b + sum(u^2)))
  sd <- sqrt((rec_eta2*(w/b + sum(u^2)))^(-1))
  
  rnorm(n = 1, mean = m, sd = sd)
}

#### Full conditional update for eta ####
update_rec_eta2 <- function(csi, beta, w, u, y, b, n_oss) {
  
  sp <- (n_oss+4)/2
  rt <- (1/2)*(((w*beta^2)/b)+sum((y-csi-beta*u)^2))
  
  rgamma(n = 1, shape = sp, rate = rt)
}

gibbs <- function(y, n_iter, n_oss, B, thin, init, prior) {
  
  ## initialize
  w_out <- numeric(n_iter)
  u_out <- matrix(NA, n_iter, n_oss)
  csi_out <- numeric(n_iter)
  beta_out <- numeric(n_iter)
  rec_eta2_out <- numeric(n_iter)
  
  beta_now <- init$beta
  rec_eta2_now <- init$rec_eta2
  csi_now <- init$csi
  
  ## Gibbs sampler
  for (i in 1:n_iter) {
    
    w_now <- update_w(beta = beta_now, 
                      rec_eta2 = rec_eta2_now, 
                      b = prior$b, 
                      d = prior$d)
    u_now <- update_u(y = y, 
                      beta = beta_now, 
                      rec_eta2 = rec_eta2_now, 
                      csi = csi_now, 
                      n_oss = n_oss)
    csi_now <- update_csi(beta = beta_now, 
                          rec_eta2 = rec_eta2_now, 
                          u = u_now, 
                          y = y, 
                          n_oss = n_oss)
    beta_now <- update_beta(csi = csi_now, 
                            rec_eta2 = rec_eta2_now, 
                            w = w_now, 
                            u = u_now, 
                            y = y, 
                            b = prior$b)
    rec_eta2_now <- update_rec_eta2(csi = csi_now, 
                                    beta = beta_now, 
                                    w = w_now, 
                                    u = u_now, 
                                    y = y, 
                                    b = prior$b, 
                                    n_oss = n_oss)
    
    w_out[i] <- w_now
    u_out[i,] <- u_now
    csi_out[i] <- csi_now
    beta_out[i] <- beta_now
    rec_eta2_out[i] <- rec_eta2_now
  }
  out <- cbind(beta = beta_out, rec_eta = rec_eta2_out, 
               csi = csi_out)[-(1:B),]
  out <- out[seq(1, nrow(out), by = thin),]
  return(out)
}

# Prior parameters
prior <- list()
prior$b <- pi^2/4
prior$d <- 0.5

# Initial values
init <- list()
init$beta <- 1 
init$rec_eta2 <- 1
init$csi <- 1



#### Tumor ####

# Data
y_1 <- read.csv("miR_182_colon.csv", sep=";")
# y_1 <- read.csv("mi_R183_colon.csv", sep=";")
# y_1 <- read.csv("mir_96_colon.csv", sep=";")
y_1 <- y_1$expression[y_1$cat == "Tum"]

set.seed(53)
R = 600000
B = 300000
post_1 <- gibbs(y = y_1, n_iter = R, B = B, 
                thin = 10, n_oss = length(y_1), 
                init = init, prior = prior)
summary(post_1)

beta_1 <- post_1[,1]
eta_1 <- (post_1[,2])^(-1/2)
csi_1 <- post_1[,3]

# Traceplot
x11()
plot(as.mcmc(post_1))

# acf
x11()
autocorr.plot(as.mcmc(post_1),lag.max = 10000)

# Plot cumulative averages
x11()
par(mfrow=c(2,2))
plot(cumsum(beta_1) / seq_along(beta_1), type = 'l', xlab = 't',
     ylab = 'X_t', main = 'Ergodic Mean plot beta_1', col="red")
plot(cumsum(eta_1) / seq_along(eta_1), type = 'l', xlab = 't',
     ylab = 'X_t', main = 'Ergodic Mean plot eta_1', col="red")
plot(cumsum(csi_1) / seq_along(csi_1), type = 'l', xlab = 't',
     ylab = 'X_t', main = 'Ergodic Mean plot csi_1', col="red")


#### Healthy ####

# Data
y_2 <- read.csv("miR_182_colon.csv", sep=";")
# y_2 <- read.csv("miR_224_colon.csv", sep=";")
# y_2 <- read.csv("mi_R183_colon.csv", sep=";")
# y_2 <- read.csv("mir_96_colon.csv", sep=";")
y_2 <- y_2$expression[y_2$cat == "Norm"]

set.seed(53)
post_2 <- gibbs(y = y_2, n_iter = R, B = B,
                thin = 10, n_oss = length(y_2), 
                init = init, prior = prior)
summary(post_2)

beta_2 <- post_2[,1]
eta_2 <- (post_2[,2])^(-1/2)
csi_2 <- post_2[,3]

# Traceplot
x11()
plot(as.mcmc(post_2))

# acf
x11()
autocorr.plot(as.mcmc(post_2),lag.max = 10000)

# Plot cumulative averages
x11()
par(mfrow=c(2,2))
plot(cumsum(beta_2) / seq_along(beta_2), type = 'l', xlab = 't',
     ylab = 'X_t', main = 'Ergodic Mean plot beta_2', col="red")
plot(cumsum(eta_2) / seq_along(eta_2), type = 'l', xlab = 't',
     ylab = 'X_t', main = 'Ergodic Mean plot eta_2', col="red")
plot(cumsum(csi_2) / seq_along(csi_2), type = 'l', xlab = 't',
     ylab = 'X_t', main = 'Ergodic Mean plot csi_2', col="red")



#### BDM ####
sigma <- function(beta, eta){
  out <- sqrt(beta^2 + eta^2)
  return(out)
}

lambda <- function(beta, eta){
  out <- beta/eta
  return(out)
}

delta <- function(beta, eta){
  out <- beta/sqrt(eta^2 + beta^2)
  return(out)
}

mu_1 <- csi_1
mu_2 <- csi_2
sigma_1 <- sigma(beta_1, eta_1)
sigma_2 <- sigma(beta_2, eta_2)
lambda_1 <- lambda(beta_1, eta_1)
lambda_2 <- lambda(beta_2, eta_2)
delta_1 <- delta(beta_1, eta_1)
delta_2 <- delta(beta_2, eta_2)

# CVs
psi_1 <- sqrt((sigma_1^2)*(1-((2/pi)*delta_1^2)))/abs(mu_1+sigma_1*delta_1*sqrt(2/pi))
psi_2 <- sqrt((sigma_2^2)*(1-((2/pi)*delta_2^2)))/abs(mu_2+sigma_2*delta_2*sqrt(2/pi))

S <- length(mu_1)

I <- sum(psi_1 < psi_2)/S

d_H <- function(Int){
  
  if(Int < 0.5){
    out <- 1 - 2*I
  }
  else {
    out <- 1 - 2*(1 - Int)
  }
}

# BDM
delta_H <- d_H(Int = I)
delta_H 

