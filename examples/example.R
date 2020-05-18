##############################################################
## Purpose: generate data from a causal spatio-temporal model
## with a one-dimensional time-invariant latent variable
##############################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(gridExtra)
library(fields)
library(foreach)
library(dplyr, warn.conflicts = FALSE)
theme_set(theme_bw())

# exponential covariance matrix
exp_cov <- function(distances, l = 1) {
  exp(-distances / l)
}

# samples spatially correlated random variables
sample_field <- function(chol_Sigma, plt = TRUE) {
  
  n <- nrow(chol_Sigma)
  z <- as.vector(t(chol_Sigma) %*% rnorm(n))
  df <- as.data.frame(coordinates)
  colnames(df) <- c("x", "y")
  df$z <- z
  
  if(plt){
    print(
      ggplot(df, aes(x, y, fill = z)) +
        geom_tile() +
        scale_fill_distiller(palette = 'Spectral') +
        coord_fixed() +
        ggtitle(as.character(substitute(cov_fun)))
    )
  }
  z
}

# samples an expontially decaying spatial trend
sample_trend <- function(coordinates=S, ell, plt = TRUE) {
  
  n <- nrow(coordinates)
  coordinates <- as.matrix(coordinates)
  dist.sq <- apply(coordinates, MARGIN = 1, function(x) sum(x^2))
  z <- exp(-dist.sq/ell)
  
  df <- as.data.frame(coordinates)
  colnames(df) <- c("x", "y")
  df$z <- z
  
  if(plt){
    print(
      ggplot(df, aes(x, y, fill = z)) +
        geom_tile() +
        scale_fill_distiller(palette = 'Spectral') +
        coord_fixed() +
        ggtitle(as.character(substitute(cov_fun)))
    )
  }
  z
}

######################  
## parameters
###################### 

N1 <- 25 # number of spatial coordinates along first
N2 <- 25 # and along second axis
N <- N1*N2 # total number of spatial locations
M <- 100 # number of time points

S <- expand.grid(x=1:N1, y = 1:N2) # spatial grid
coordinates <- as.matrix(S)
distances <- fields::rdist(coordinates) # pairwise distances between points
Sigma_exp <- exp_cov(distances, l = 2) + 1e-6 # covariance matrix
diag(Sigma_exp) <- diag(Sigma_exp) + 1e-6 # avoiding numerical underflow
chol_Sigma_exp <- chol(Sigma_exp) # cholesky decomposition

muH1 <- 0 # mean of H1
sH1 <- 1 # sd of H1
muH2 <- 1 # mean of H2
sH2 <- 1 # sd of H2
rho <- .5 # product mean E(H1*H2) 


muS_X <- sample_trend(S, 1e3, TRUE) # spatial trend in X
muT_X <- 1 + .5*sin(2*pi*(1:M)/M) # temporal trend in X

set.seed(1)
X <- Y <- H1 <- H2 <- matrix(NA,M,N)
h1 <- sample_field(chol_Sigma_exp,plt=TRUE)
h2 <- muH2 + rho*h1 + sqrt(1-rho^2)*sample_field(chol_Sigma_exp,plt=TRUE)
for(t in 1:M){
  H1[t,] <- h1
  H2[t,] <- h2
  X[t,] <- muS_X + muT_X[t]*H1[t,]*H2[t,]/5 + .5*sample_field(chol_Sigma_exp,plt=FALSE)
  Y[t,] <- X[t,]*(3/2+H1[t,]*H2[t,]) + H1[t,]^2 + abs(H2[t,])*sample_field(chol_Sigma_exp,plt=FALSE)
}
############################################################
## GROUND TRUTH: 
## f_AVE(X->Y)(x) = (3/2 + E(H1*H2))*x + E(H1^2) = (3/2 + rho)*x + sH1 = 2*x + 1
b1 <- rho + 3/2
b0 <- sH1 
############################################################

######################  
## estimation
###################### 
# data ordering: (s1,t1), ..., (s1,tM), ..., (sN,t1), ..., (sN,tM)
# OBS: important for estimation and for resampling procedure
X <- c(X) 
Y <- c(Y)
H1 <- c(H1)
H2 <- c(H2)


# estimated regression coefficients in each separate location
bhat <- function(X,Y){
  foreach(s = 1:N, .combine = "rbind") %do% {
    #print(paste(i, "out of ", N1*N2))
    Xs <- X[(M*(s-1)+1):((M*(s-1)+M))]
    Ys <- Y[(M*(s-1)+1):((M*(s-1)+M))]
    coefficients(lm(Ys~Xs))
  }
}

# estimated causal coefficients
b.causal <- function(X,Y) colMeans(bhat(X,Y))
bhat.data <- bhat(X,Y)


# estimate of causal coefficients (true values = (1,2))
(b.causal.data <- b.causal(X,Y)) # (0.88, 1.93)

######################  
## plotting
###################### 

linedata <- data.frame(intercept = bhat.data[,1], 
                       slope = bhat.data[,2])


# red: nonparametric regression
# green: avg causal effect (unknown)
# thin blue: linear regressions within each locatino
# thick blue: estimate of the avg causal effect
p <- ggplot(df, aes(X,Y)) + geom_point(alpha=.7, aes(col=H1*H2)) + 
  scale_color_gradient(name = expression(paste(bar(H)," * ",tilde(H))), low="black", high="gray") + 
  geom_abline(data=linedata, aes(slope=slope, intercept=intercept), col = "blue", alpha = .15, size = .5) + 
  geom_abline(intercept = b.causal.data[1], slope = b.causal.data[2], col = "blue", size=1.5) + 
  theme(legend.position = c(.2,.7), 
        text = element_text(size=13)) + 
  geom_smooth(col = "red", se = FALSE, size = 1.5) + 
  geom_abline(intercept = b0[1], slope = b0[1], col = "green", size = 1.5, lty=2)
p


########################
## resampling test
########################

# resampling procedure
res <- function(Y){
  z <- sample(1:M, size = M, replace=FALSE) # temporal permutation
  res.ind <- M*rep(0:(N-1), each=M) + rep(z, N) # same temporal permutation in each spatial location
  Y[res.ind]
}

# test statistic
ts <- function(X,Y) b.causal(X,Y)[2]

# two-sided p-value
B <- 99 # number of permutations
pval <- function(X,Y,resampler=res,test.stat=ts,plt=TRUE){
  t <- test.stat(X,Y)
  t.res <- foreach(b=1:B, .combine = "c") %do% {
    print(b)
    set.seed(b)
    Y.res <- resampler(Y)
    test.stat(X,Y.res)
  }
  pv <- min(1,2*min((1+sum(t.res <= t))/(1+B), (1+sum(t.res >= t))/(1+B)))
  if(plt){
    dev.off()
    hist(t.res, breaks = 20, xlim = c(-1,1)*1.2*abs(t))
    abline(v=t, col="red")
  }
  pv
}

(pv <- pval(X,Y,res,ts)) # pv = 0.02 (smallest attainable p-value)