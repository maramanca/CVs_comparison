#### NB example ####

########################## DATA Population 1 ##########################################################
### Laxminarayan et al. "Epidemiology and transmission dynamics of COVID-19 in two Indian states".  ### 
### In: Science(2020). doi:https://www.science.org/doi/10.1126/science.abd7672                      ###
#######################################################################################################

########################## DATA Population 2 ##########################################################                                                                               
### Adam et al. "Clustering and superspreading potential of SARS-C-oV-2 in-fections in Hong Kong".  ###
### In: Nature Medicine(2020). doi:https://doi.org/10.1038/s41591-020-1092-0                        ###
#######################################################################################################

library(mvtnorm)

#### Reading the data of the two populations from .csv files and setting the population parameters ####

Pop1 <- read.csv("dati_india.csv")
Pop2 <- read.csv("dati_hongkong.csv")

x1 <- Pop1[, 2]
x2 <- Pop2[, 2]

n1 <- length(x1)
m1 <- mean(x1)

n2 <- length(x2)
m2 <- mean(x2)

# Function involved in the expression of the log-posterior distribution
h <- function(x, theta){
  
  out <- rep(NA, length(x))
  for (i in 1:length(x)) {
    out[i] <- factorial(x[i] + theta[1] - 1)
  }
  
  prod <- sum(log(out))
  
  return(prod)
}

# Negative Binomial log-posterior distribution with parameter theta = (alpha, beta)
log_posterior_NB <- function(theta, x, n, m){
  
  pd.1 <- h(x, theta)
  pd.2 <- (-n)*log(factorial(theta[1] - 1))
  pd.3 <- n*theta[1]*log(theta[2]) + log(theta[1])
  pd.4 <- -n*(theta[1] + m)*log(theta[2] + 1)
  pd.5 <- (1/2)*log(theta[1]*psigamma(theta[1], 1) - 1)
  
  out <- pd.1 + pd.2 + pd.3 + pd.4 + pd.5
  
  return(out)
}



#### Metropolis-Hastings algorithm for the Negative Binomial ####

set.seed(5)
MH_NB <- function(S, B, theta0, x, n, m, tau){
  
  R <- B + S
  mat <- matrix(NA, nrow = R, ncol = length(theta0))
  acc <- rep(0, R)
  
  theta_chain <- theta0
  
  for (i in 1:R) {
    # Sampling from the proposal 
    l_theta_star <- rmvnorm(1, mean = log(theta_chain), sigma = tau)
    theta_star <- exp(l_theta_star)
    
    # Computation of r
    l1 <- log_posterior_NB(theta = theta_star, x = x, n = n, m = m)
    l2 <- log_posterior_NB(theta = theta_chain, x = x, n = n, m = m)
    r <- min(1, exp(l1 - l2))
    
    # Drawing U from the uniform distribution in (0,1)
    u <- runif(1)
    
    # Accepting or rejecting theta_star
    if(r >= u) {
      theta_chain <- theta_star
      acc[i] <- 1
    }
    mat[i, ] <- theta_chain
  }
  
  # Output from discarding burn-in values
  list(values = mat[(B + 1):R,], acc_rate=sum(acc[(B + 1):R])/S)
}



#### Application of the sampling algorithm to Population 1 ####

# theta_0 parameter initialization
maximization_1 <- optim(function(r) -log_posterior_NB(r, x = x1, n = n1, m = m1),
                        hessian = T, method = "L-BFGS-B", 
                        par = c(2,2),
                        lower = c(0.001,0.0001), upper = c(100,100))
theta0_1 <- maximization_1$par
theta0_1 

# tau covariance matrix for the proposal distribution
tau_1 <- solve(maximization_1$hessian)
tau_1

# Sampling procedure
S <- 100000 # number of samples
B <- 30000  # burn-in samples

set.seed(1)
posterior_sample_1 <- MH_NB(S = S, B = B, theta0 = theta0_1, 
                            x = x1, n = n1, 
                            m = m1, tau = tau_1)

alpha_1 <- posterior_sample_1$values[, 1]
beta_1 <- posterior_sample_1$values[, 2]

### Convergence diagnostics ###
# Samples visualization, alpha_1 vs beta_1
plot(posterior_sample_1$values, xlab = "alpha_1", ylab = "beta_1")

# Acceptance rate
posterior_sample_1$acc_rate

# Traceplots
par(mfrow = c(1,2))
plot(alpha_1, type = "l")
plot(beta_1, type = "l")

# Histograms
means <- c(mean(alpha_1), mean(beta_1))
means
medians <- c(median(alpha_1), median(beta_1))
medians

par(mfrow = c(1,2))
{hist(alpha_1, 3000, freq = F, xlab = "alpha_1", main = "a) Histogram of alpha_1")
  abline(v = means[1], col = "red", lwd = 2)
  abline(v = medians[1], col = "green", lwd = 2, lty = 3)}
{hist(beta_1, 3000, freq = F, xlab = "beta_1", main = "b) Histogram of beta_1") 
  abline(v = means[2], col = "red", lwd = 2)
  abline(v = medians[2], col = "green", lwd = 2, lty = 3)}

# ACF
par(mfrow = c(1,2))
acf(alpha_1, main = "ACF for alpha_1", lag.max = 10000)
acf(beta_1, main = "ACF for beta_1", lag.max = 10000)



#### Application of the sampling algorithm to Population 2 ####

# theta_0 parameter initialization
maximization_2 <- optim(function(r) -log_posterior_NB(r, x = x2, n = n2, m = m2),
                        hessian = T, method = "L-BFGS-B", 
                        par = c(2,2),
                        lower = c(0.001, 0.0001), upper = c(100, 100))
theta0_2 <- maximization_2$par
theta0_2 

# tau covariance matrix for the proposal distribution
tau_2 <- 3*solve(maximization_2$hessian)
tau_2

# Sampling procedure
S <- 100000 # number of samples
B <- 30000  # burn-in samples

set.seed(1)
posterior_sample_2 <- MH_NB(S = S, B = B, theta0 = theta0_2, 
                            x = x2, n = n2, 
                            m = m2, tau = tau_2)

alpha_2 <- posterior_sample_2$values[, 1]
beta_2 <- posterior_sample_2$values[, 2]

### Convergence diagnostics ###
# Samples visualization, alpha_2 vs beta_2
plot(posterior_sample_2$values, xlab = "alpha_2", ylab = "beta_2")

# Acceptance rate
posterior_sample_2$acc_rate

# Traceplots
par(mfrow = c(1,2))
plot(alpha_2, type = "l")
plot(beta_2, type = "l")

# Histograms
means <- c(mean(alpha_2), mean(beta_2))
means
medians <- c(median(alpha_2), median(beta_2))
medians

par(mfrow = c(1,2))
{hist(alpha_1, 3000, freq = F, xlab = "alpha_1", main = "a) Histogram of alpha_2")
  abline(v = means[1], col = "red", lwd = 2)
  abline(v = medians[1], col = "green", lwd = 2, lty = 3)}
{hist(beta_1, 3000, freq = F, xlab = "beta_1", main = "b) Histogram of beta_2") 
  abline(v = means[2], col = "red", lwd = 2)
  abline(v = medians[2], col="green", lwd = 2, lty = 3)}

# ACF
par(mfrow = c(1,2))
acf(alpha_2, main = "ACF for alpha_1" , lag.max = 10000)
acf(beta_2, main = "ACF for beta_1", lag.max = 10000)



#### BDM computation and hypothesis testing H: psi_1 = psi_2 ####

psi_1 <- sqrt((beta_1 + 1)/alpha_1) 
psi_2 <- sqrt((beta_2 + 1)/alpha_2)

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


