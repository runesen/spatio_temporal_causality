setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(gridExtra)
library(fields)
library(foreach)
library(dplyr, warn.conflicts = FALSE)

theme_set(theme_bw())


exp_cov <- function(distances, l = 2) { # CHANGE back to 2
  exp(-distances / l)
}


sample_field <- function(chol_Sigma, df.coordinates, dimensions = "space", plt = TRUE) {
  
  n <- nrow(chol_Sigma)
  z <- as.vector(t(chol_Sigma) %*% rnorm(n))
  
  if(plt){
    if(dimensions == "space"){
      colnames(df.coordinates) <- c("s1", "s2")
      df.coordinates$value <- z
      p <- ggplot(df.coordinates, aes(s1, s2, fill = value)) +
        geom_tile() +
        scale_fill_distiller(palette = 'Spectral') +
        coord_fixed()
    }
    if(dimensions == "time"){
      colnames(df.coordinates) <- "t"
      df.coordinates$value <- z
      p <- ggplot(df.coordinates, aes(t, value)) + geom_line()
    }
    if(dimensions == "spacetime"){
      colnames(df.coordinates) <- c("s1", "s2", "t")
      df.coordinates$value <- z
      p <- ggplot(df.coordinates, aes(s1, s2, fill = value)) +
        geom_tile() + 
        # facet_wrap(t ~ ., ncol = 5) + 
        facet_grid(. ~t) + 
        scale_fill_distiller(palette = 'Spectral') +
        coord_fixed()
    }
    print(p)
  }
  z
}


set.seed(1)
n.lon <- 60   
# n.lon <- 6 # CHANGE 
n.lat <- 25
n.space <- n.lon*n.lat
n.time <- 50
n.rep <- 100
# n.rep <- 50 # CHANGE
sigma.seq <- sqrt(seq(0,.2,length.out = 5))
n.sigma <- length(sigma.seq) 

df.space <- expand.grid(x=1:n.lon, y = 1:n.lat) # spatial grid
df.time <- data.frame(t = 1:n.time) # temporal grid
df.spacetime <- expand.grid(s1=1:n.lon, s2 = 1:n.lat, t = 1:n.time) # spatio-temporal grid
distances.space <- fields::rdist(df.space)
distances.time <- fields::rdist(df.time)

# Exponential covariance matrix with l = 2
Sigma_spatial_error <- exp_cov(distances.space,l=2) 
diag(Sigma_spatial_error) <- diag(Sigma_spatial_error) + 1e-6
chol_Sigma_spatial_error <- chol(Sigma_spatial_error)
spatial_error <- sample_field(chol_Sigma_spatial_error, df.space, dimensions = "space", plt=TRUE)

spatial_trend_X <- rep(exp(-apply(df.space, MARGIN = 1, function(x) sum(x^2)) / 1000), n.time)
temporal_trend_X <- rep(.2 + .1*sin(2*pi*(1:n.time)/n.time), each = n.space)


# ground truth: E[Y; do(X=x)] = (1.5 + E[Hbar * Htilde]) * x + E[Hbar^2] = 2*x + 1

df <- foreach(sigma = sigma.seq, .combine = "rbind") %do% {
  a <- sqrt(1-sigma^2)
  b <- 1/(2*a)
  c <- sqrt(1-b^2)
  
  foreach(rep = 1:n.rep, .combine = "rbind") %do% {
    print(paste("sigma = ", sigma, ", rep = ", rep))
    
    zeta <- rep(sample_field(chol_Sigma_spatial_error, df.space, dimensions = "space", plt=FALSE), n.time) # zeta 
    psi <- rep(sample_field(chol_Sigma_spatial_error, df.space, dimensions = "space", plt=FALSE), n.time) # psi
    chi <- rnorm(n.time*n.space) # chi
    xi <- c(sapply(1:n.time, function(t) sample_field(chol_Sigma_spatial_error, df.space, dimensions = "space", plt=FALSE))) # xi
    epsilon <- c(sapply(1:n.time, function(t) sample_field(chol_Sigma_spatial_error, df.space, dimensions = "space", plt=FALSE))) # epsilon
    
    Hbar <- a * zeta + sigma * chi # contains the time-varying component chi
    Htilde <- 1 + b * zeta + c * psi
    X <- spatial_trend_X + temporal_trend_X * Hbar * Htilde + 0.5 * xi
    Y <- (1.5 + Hbar * Htilde) * X + Hbar^2 + abs(Htilde) * epsilon
    
    tmp <- df.spacetime
    tmp$rep <- rep
    tmp$sigma <- sigma
    tmp$Hbar <- Hbar 
    tmp$Htilde <- Htilde
    tmp$X <- X
    tmp$Y <- Y
    tmp
  }
}


# all settings for which we require parameter estimates of b0 and b1
df.coef <- expand.grid(sigma = sigma.seq, rep = 1:n.rep)
n.settings <- nrow(df.coef)
df.coef <- df.coef[rep(1:n.settings,3),]
df.coef$method <- rep(c("marginal", "adj-for-time-inv-H", "adj-for-all-H"), each=n.settings)
df.coef$b0 <- df.coef$b1 <- NA


#######################
## marginal regression
#######################
for(i in 1:n.settings){
  print(i)
  df.i <- subset(df, sigma == df.coef$sigma[i] & rep == df.coef$rep[i])
  bhat <- coef(lm(Y~X, data=df.i))
  df.coef$b0[i] <- bhat[1]
  df.coef$b1[i] <- bhat[2]
}

#######################
## adjusting for H
#######################

# sort s.t. the most inner loop is time
df <- df[order(df$sigma,df$rep,df$s1,df$s2,df$t),] 
hbar <- df$Hbar
htilde <- df$Htilde
x <- df$X
y <- df$Y
  
# for each setting, at each spatial location, do: 
# * a marginal regression of Ys onto Xs, AND
# * a regression of Ys onto (Xs, Hbars) & averaging out Hbar
coef.per.location <- foreach(i = 1:(n.rep*n.sigma), .combine = "rbind") %do% {
  index.dataset <- n.time*n.space*(i-1) # last index before current dataset
  hbar.i <- hbar[(index.dataset+1):(index.dataset + n.time*n.space)] # realization of Hbar for current setting and repetition
  foreach(k = 1:n.space, .combine = "rbind") %do% {
    print(paste("setting ", i, "out of ", n.rep*n.sigma, ", location ", k, " out of ", n.space))
    index.location <- index.dataset + n.time*(k-1) # last index before currect location
    foreach(j = (n.time/10), .combine = "rbind") %do% {
      # data from setting i in spatial location k and time points 1, ..., 10*j
      x.jk <- x[(index.location+1):(index.location+j*10)]
      y.jk <- y[(index.location+1):(index.location+j*10)]
      mod.adj.time.inv.H <- lm(y.jk ~ x.jk) # regression of Y onto X

      hbar.jk <- hbar[(index.location+1):(index.location+j*10)]
      mod.adj.all.H <- lm(y.jk ~ -1 + x.jk + I(hbar.jk^2) + x.jk:hbar.jk) # regression of Y onto (X,H)
      
      c(coef(mod.adj.time.inv.H), # estimates of b0 and b1 without adjusting for Hbar
        coef(mod.adj.all.H)[2]*mean(hbar.jk^2), # estimate of b0 when adjusting for Hbar (and then averaging out Hbar)
        coef(mod.adj.all.H)[1] + coef(mod.adj.all.H)[3]*mean(hbar.jk))  # estimate of b1 when adjusting for Hbar (and then averaging out Hbar)
    }
  }
}



df.sub <- subset(df, t == n.time)
df.sub <- df.sub[rep(1:nrow(df.sub),2),] # duplicating
df.sub$b0 <- c(coef.per.location[,1], coef.per.location[,3])
df.sub$b1 <- c(coef.per.location[,2], coef.per.location[,4])
df.sub$method <- factor(rep(c("adj-for-time-inv-H", "adj-for-all-H"), each = nrow(df.sub)/2), levels = c("adj-for-time-inv-H", "adj-for-all-H"))
df.sub <- df.sub[order(df.sub$method, df.sub$sigma, df.sub$rep, df.sub$t, df.sub$s1, df.sub$s2),] # inner loop is latitude, then longitude, then time
b0hat <- df.sub$b0
b1hat <- df.sub$b1


bhat <- foreach(i = 1:(2*n.sigma*n.rep), .combine = "rbind") %do% { # loop over methods, settings and replications
  print(paste(i, "out of ", 2*n.sigma*n.rep))
  foreach(j = (n.lon/10), .combine = "rbind") %do% { # for increasing spatial domain, compute spatial averages
    index <- (n.space*(i-1)+1):(n.space*(i-1)+n.lat*j*10)
    c(mean(b0hat[index],na.rm=1), mean(b1hat[index],na.rm=1))
  }
}



df.coef$method <- factor(df.coef$method, levels = c("adj-for-time-inv-H", "adj-for-all-H", "marginal"))
df.coef <- df.coef[order(df.coef$method, df.coef$sigma, df.coef$rep),]
df.coef$b0[1:(2*n.settings)] <- bhat[,1]
df.coef$b1[1:(2*n.settings)] <- bhat[,2]
df.coef$err <- sqrt((df.coef$b0-1)^2 + (df.coef$b1-2)^2)
df.coef$var <- df.coef$sigma^2
write.table(df.coef, "time-var-conf.txt", quote = FALSE)



