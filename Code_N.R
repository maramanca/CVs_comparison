#### N example ####

library(NormalGamma)

# Posterior distribution
rnormgamma <- function(n, eta, nu, alpha, beta) {
  
  if (length(n) > 1) 
  n <- length(n)
  phi <- rgamma(n, alpha, beta)
  mu <- rnorm(n, eta, 1/sqrt(nu*phi))
  
  data.frame(mu = mu, phi = phi)
}

# Function that gives the hyperparameters of the posterior distribution
hyperparam <- function(n, mean, s_square){
  eta <- mean          
  nu <- n
  alpha <- (n - 1)/2
  beta <- (n*s_square)/2
  list(eta = eta, nu = nu, alpha = alpha, beta = beta)
} 

# Function that computes the BDM 
d_H <- function(Int){
  if(Int < 0.5){
    out <- 1 - 2*Int
  }
  else {
    out <- 1 - 2*(1 - Int)
  }
}

#### Weight ####
# [1]
n_1 <- 140
x1_mean <- 67.22
s_square1 <- 8.46^2
# [2]
n_2 <- 140
x2_mean <- 53.71
s_square2 <- 7.59^2

# Hyperparameters of the two populations
h1 <- hyperparam(n = n_1, mean = x1_mean, s_square = s_square1)
h2 <- hyperparam(n = n_2, mean = x2_mean, s_square = s_square2)

# Sampling of the parameter vector theta = (phi,mu) for the two populations
set.seed(1)
m <- 10000
phi_mu_1 <- rnormgamma(n = m, 
                       eta = h1$eta, nu = h1$nu, 
                       alpha = h1$alpha, beta = h1$beta)
phi_mu_2 <- rnormgamma(n = m, 
                       eta = h2$eta, nu = h2$nu, 
                       alpha = h2$alpha, beta = h2$beta)

### BDM computation and hypothesis testing H: psi_1 = psi_2 ###
psi_1 <- 1/(phi_mu_1$mu*(phi_mu_1$phi)^(0.5))
psi_2 <- 1/(phi_mu_2$mu*(phi_mu_2$phi)^(0.5))

I <- sum(psi_1<psi_2)/m

delta_H <- d_H(Int = I)
delta_H
# 0.812 do not reject H0

#### Cephalic breadth ####
# [1]
n_1 <- 141
x1_mean <- 15.10
s_square1 <- 0.64^2
# [2]
n_2 <- 172
x2_mean <- 14.53
s_square2 <- 0.58^2

# Hyperparameters of the two populations
h1 <- hyperparam(n = n_1, mean = x1_mean, s_square = s_square1)
h2 <- hyperparam(n = n_2, mean = x2_mean, s_square = s_square2)

# Sampling of the parameter vector theta = (phi,mu) for the two populations
set.seed(1)
m <- 10000
phi_mu_1 <- rnormgamma(n = m, 
                       eta = h1$eta, nu = h1$nu, 
                       alpha = h1$alpha, beta = h1$beta)
phi_mu_2 <- rnormgamma(n = m, 
                       eta = h2$eta, nu = h2$nu, 
                       alpha = h2$alpha, beta = h2$beta)

### BDM computation and hypothesis testing H: psi_1 = psi_2 ###
psi_1 <- 1/(phi_mu_1$mu*(phi_mu_1$phi)^(0.5))
psi_2 <- 1/(phi_mu_2$mu*(phi_mu_2$phi)^(0.5))

I <- sum(psi_1 < psi_2)/m

# BDM
delta_H <- d_H(Int = I)
delta_H
# 0.550 do not reject H0

#### Elbow ####
# [1]
n_1 <- 103
x1_mean <- 7.02
s_square1 <- 0.39^2
# [2]
n_2 <- 117
x2_mean <- 6.01
s_square2 <- 0.35^2

# Hyperparameters of the two populations
h1 <- hyperparam(n = n_1, mean = x1_mean, s_square = s_square1)
h2 <- hyperparam(n = n_2, mean = x2_mean, s_square = s_square2)

# Sampling of the parameter vector theta = (phi,mu) for the two populations
set.seed(1)
m <- 10000
phi_mu_1 <- rnormgamma(n = m, 
                       eta = h1$eta, nu = h1$nu, 
                       alpha = h1$alpha, beta = h1$beta)
phi_mu_2 <- rnormgamma(n = m, 
                       eta = h2$eta, nu = h2$nu, 
                       alpha = h2$alpha, beta = h2$beta)

### BDM computation and hypothesis testing H: psi_1 = psi_2 ###
psi_1 <- 1/(phi_mu_1$mu*(phi_mu_1$phi)^(0.5))
psi_2 <- 1/(phi_mu_2$mu*(phi_mu_2$phi)^(0.5))

I <- sum(psi_1<psi_2)/m

delta_H  <-  d_H(Int = I)
delta_H
# 0.355 do not reject H0

#### Midarm relaxed ####
# [1]
n_1 <- 139
x1_mean <- 26.91
s_square1 <- 2.60^2
# [2]
n_2 <- 134
x2_mean <- 23.47
s_square2 <- 2.01^2

# Hyperparameters of the two populations
h1 <- hyperparam(n = n_1, mean = x1_mean, s_square = s_square1)
h2 <- hyperparam(n = n_2, mean = x2_mean, s_square = s_square2)

# Sampling of the parameter vector theta = (phi,mu) for the two populations
set.seed(1)
m <- 10000
phi_mu_1 <- rnormgamma(n = m, 
                       eta = h1$eta, nu = h1$nu, 
                       alpha = h1$alpha, beta = h1$beta)
phi_mu_2 <- rnormgamma(n = m, 
                       eta = h2$eta, nu = h2$nu, 
                       alpha = h2$alpha, beta = h2$beta)

### BDM computation and hypothesis testing H: psi_1 = psi_2 ###
psi_1 <- 1/(phi_mu_1$mu*(phi_mu_1$phi)^(0.5))
psi_2 <- 1/(phi_mu_2$mu*(phi_mu_2$phi)^(0.5))

I <- sum(psi_1 < psi_2)/m

# BDM
delta_H <- d_H(Int = I)
delta_H
# 0.831 do not reject H0

#### Midarm tensed ####
# [1]
n_1 <- 139
x1_mean <- 30.83
s_square1 <- 2.74^2
# [2]
n_2 <- 133
x2_mean <- 25.28
s_square2 <- 2.15^2

# Hyperparameters of the two populations
h1 <- hyperparam(n = n_1, mean = x1_mean, s_square = s_square1)
h2 <- hyperparam(n = n_2, mean = x2_mean, s_square = s_square2)

# Sampling of the parameter vector theta = (phi,mu) for the two populations
set.seed(1)
m <- 10000
phi_mu_1 <- rnormgamma(n = m, 
                       eta = h1$eta, nu = h1$nu, 
                       alpha = h1$alpha, beta = h1$beta)
phi_mu_2 <- rnormgamma(n = m, 
                       eta = h2$eta, nu = h2$nu, 
                       alpha = h2$alpha, beta = h2$beta)

### BDM computation and hypothesis testing H: psi_1 = psi_2 ###
psi_1 <- 1/(phi_mu_1$mu*(phi_mu_1$phi)^(0.5))
psi_2 <- 1/(phi_mu_2$mu*(phi_mu_2$phi)^(0.5))

I <- sum(psi_1<psi_2)/m

delta_H <- d_H(Int = I)
delta_H
# 0.388 do not reject H0

#### Biceps ####
# [1]
n_1 <- 137
x1_mean <- 4.10
s_square1 <- 1.79^2
# [2]
n_2 <- 133
x2_mean <- 6.08
s_square2 <- 2.51^2

# Hyperparameters of the two populations
h1 <- hyperparam(n = n_2, mean = x1_mean, s_square = s_square1)
h2 <- hyperparam(n = n_2, mean = x2_mean, s_square = s_square2)

# Sampling of the parameter vector theta = (phi,mu) for the two populations
set.seed(1)
m <- 10000
phi_mu_1 <- rnormgamma(n = m, 
                       eta = h1$eta, nu = h1$nu, 
                       alpha = h1$alpha, beta = h1$beta)
phi_mu_2 <- rnormgamma(n = m, 
                       eta = h2$eta, nu = h2$nu, 
                       alpha = h2$alpha, beta = h2$beta)

### BDM computation and hypothesis testing H: psi_1 = psi_2 ###
psi_1 <- 1/(phi_mu_1$mu*(phi_mu_1$phi)^(0.5))
psi_2 <- 1/(phi_mu_2$mu*(phi_mu_2$phi)^(0.5))

I <- sum(psi_1 < psi_2)/m

# BDM
delta_H <- d_H(Int = I)
delta_H
# 0.4202 do not reject H0


#### Triceps ####
# [1]
n_1 <- 140
x1_mean <- 7.76
s_square1 <- 3.76^2
# [2]
n_2 <- 140
x2_mean <- 13.33
s_square2 <- 4.78^2

# Hyperparameters of the two populations
h1 <- hyperparam(n = n_1, mean = x1_mean, s_square = s_square1)
h2 <- hyperparam(n = n_2, mean = x2_mean, s_square = s_square2)

# Sampling of the parameter vector theta = (phi,mu) for the two populations
set.seed(1)
m <- 10000
phi_mu_1 <- rnormgamma(n = m, 
                       eta = h1$eta, nu = h1$nu, 
                       alpha = h1$alpha, beta = h1$beta)
phi_mu_2 <- rnormgamma(n = m, 
                       eta = h2$eta, nu = h2$nu, 
                       alpha = h2$alpha, beta = h2$beta)

### BDM computation and hypothesis testing H: psi_1 = psi_2 ###
psi_1 <- 1/(phi_mu_1$mu*(phi_mu_1$phi)^(0.5))
psi_2 <- 1/(phi_mu_2$mu*(phi_mu_2$phi)^(0.5))

I <- sum(psi_1<psi_2)/m

delta_H <- d_H(Int = I)
delta_H
# 0.9962 reject H0


#### Subscapular ####
# [1]
n_1 <- 137
x1_mean <- 10.34
s_square1 <- 3.78^2
# [2]
n_2 <- 140
x2_mean <- 12.71
s_square2 <- 4.53^2

# Hyperparameters of the two populations
h1 <- hyperparam(n = n_1, mean = x1_mean, s_square = s_square1)
h2 <- hyperparam(n = n_2, mean = x2_mean, s_square = s_square2)

# Sampling of the parameter vector theta = (phi,mu) for the two populations
set.seed(1)
m <- 10000
phi_mu_1 <- rnormgamma(n = m, 
                       eta = h1$eta, nu = h1$nu, 
                       alpha = h1$alpha, beta = h1$beta)
phi_mu_2 <- rnormgamma(n = m, 
                       eta = h2$eta, nu = h2$nu, 
                       alpha = h2$alpha, beta = h2$beta)

### BDM computation and hypothesis testing H: psi_1 = psi_2 ###
psi_1 <- 1/(phi_mu_1$mu*(phi_mu_1$phi)^(0.5))
psi_2 <- 1/(phi_mu_2$mu*(phi_mu_2$phi)^(0.5))

I <- sum(psi_1 < psi_2)/m

# BDM
delta_H <- d_H(Int = I)
delta_H
# 0.213 do not reject H0


#### Suprailiac ####
# [1]
n_1 <- 140
x1_mean <- 9.23
s_square1 <- 4.34^2
# [2]
n_2 <- 140
x2_mean <- 10.21
s_square2 <- 4.48^2

# Hyperparameters of the two populations
h1 <- hyperparam(n = n_1, mean = x1_mean, s_square = s_square1)
h2 <- hyperparam(n = n_2, mean = x2_mean, s_square = s_square2)

# Sampling of the parameter vector theta = (phi,mu) for the two populations
set.seed(1)
m <- 10000
phi_mu_1 <- rnormgamma(n = m, 
                       eta = h1$eta, nu = h1$nu, 
                       alpha = h1$alpha, beta = h1$beta)
phi_mu_2 <- rnormgamma(n = m, 
                       eta = h2$eta, nu = h2$nu, 
                       alpha = h2$alpha, beta = h2$beta)

### BDM computation and hypothesis testing H: psi_1 = psi_2 ###
psi_1 <- 1/(phi_mu_1$mu*(phi_mu_1$phi)^(0.5))
psi_2 <- 1/(phi_mu_2$mu*(phi_mu_2$phi)^(0.5))

I <- sum(psi_1 < psi_2)/m

# BDM
delta_H <- d_H(Int = I)
delta_H
# 0.5072 do not reject H0


#### Abdominal ####
# [1]
n_1 <- 97
x1_mean <- 12.15
s_square1 <- 6.52^2
# [2]
n_2 <- 111
x2_mean <- 12.77
s_square2 <- 5.74^2

# Hyperparameters of the two populations
h1 <- hyperparam(n = n_1, mean = x1_mean, s_square = s_square1)
h2 <- hyperparam(n = n_2, mean = x2_mean, s_square = s_square2)

# Sampling of the parameter vector theta = (phi,mu) for the two populations
set.seed(1)
m <- 10000
phi_mu_1 <- rnormgamma(n = m, 
                       eta = h1$eta, nu = h1$nu, 
                       alpha = h1$alpha, beta = h1$beta)
phi_mu_2 <- rnormgamma(n = m, 
                       eta = h2$eta, nu = h2$nu, 
                       alpha = h2$alpha, beta = h2$beta)

### BDM computation and hypothesis testing H: psi_1 = psi_2 ###
psi_1 <- 1/(phi_mu_1$mu*(phi_mu_1$phi)^(0.5))
psi_2 <- 1/(phi_mu_2$mu*(phi_mu_2$phi)^(0.5))

I <- sum(psi_1 < psi_2)/m

# BDM
delta_H <- d_H(Int = I)
delta_H
# 0.848 do not reject H0
