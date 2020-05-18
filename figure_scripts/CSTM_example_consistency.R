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


sq_exp_cov <- function(distances, l = 1) {
  exp(-distances ^ 2 / l ^ 2)
}



exp_cov <- function(distances, l = 1) {
  exp(-distances / l)
}

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


############## samples ##############
set.seed(1)
N1 <- 100
N2 <- 25
N <- N1*N2
M <- 100
B <- 100

S <- expand.grid(x=1:N1, y = 1:N2) # spatial grid
coordinates <- as.matrix(S)
distances <- fields::rdist(coordinates)
# Sigma_sq_exp <- sq_exp_cov(distances)
# diag(Sigma_sq_exp) <- diag(Sigma_sq_exp) + 1e-6
# chol_Sigma_sq_exp <- chol(Sigma_sq_exp)
Sigma_exp <- exp_cov(distances, l = 2) + 1e-6
diag(Sigma_exp) <- diag(Sigma_exp) + 1e-6
chol_Sigma_exp <- chol(Sigma_exp)
H1 <- sample_field(chol_Sigma_exp,plt=TRUE);mean(H1)

muH1 <- 0
sH1 <- 1
muH2 <- 1
sH2 <- 1
rho <- .5 # cov(H1, H2) if muH1 = muH2 = 0 and sH1 = 1

## ground truth: 
# AVG causal effect: f(x) = E[H1*H2]*x + E[H1^2] = rho*x + sH1 
b0 <- sH1
b1 <- rho + 3/2

S <- expand.grid(x=1:N1, y = 1:N2) # spatial grid

muS_X <- sample_trend(S, 1e3, TRUE) # spatial trend in X
muT_X <- 1 + .5*sin(2*pi*(1:M)/M);plot(muT_X)


H1 <- sample_field(chol_Sigma_exp,plt=TRUE)
H2 <- muH2 + rho*H1 + sqrt(1-rho^2)*sample_field(chol_Sigma_exp,plt=TRUE)

if(0){
  
  s1vec <- expand.grid(s1=1:N1, s2=1:N2)$s1
  b1hat <- sapply(1:N1, function(ss1) mean(H1[s1vec<=ss1]*H2[s1vec<=ss1]))
  plot(b1hat)
  abline(h=.5,lwd=3,col="blue")
  lines(1:N1, .5-1.96/sqrt(N2*(1:N1)),col="red",lwd=3)
  lines(1:N1, .5+1.96/sqrt(N2*(1:N1)),col="red",lwd=3)
  
  b0hat <- sapply(1:N1, function(ss1) mean(H1[s1vec<=ss1]^2))
  plot(b0hat)
  abline(h=1,lwd=3,col="blue")
  lines(1:N1, 1-1.96/sqrt(N2*(1:N1)),col="red",lwd=3)
  lines(1:N1, 1+1.96/sqrt(N2*(1:N1)),col="red",lwd=3)
  
  muRhat <- sapply(1:N, function(i) mean(R1[1:i]))
  plot(muRhat)
  abline(h=0,lwd=3,col="blue")
  lines(1:N, -1.96/(sqrt(1:N)),lwd=3,col="red")
  lines(1:N, +1.96/(sqrt(1:N)),lwd=3,col="red")
}

x.out <- y.out <- h1.out <- h2.out <- c()
for(i in 1:B){
  set.seed(i)
  print(paste0("b = ", i))
  X <- Y <- H1 <- H2 <- list()
  h1 <- sample_field(chol_Sigma_exp,plt=FALSE)
  h2 <- muH2 + rho*h1 + sqrt(1-rho^2)*sample_field(chol_Sigma_exp,plt=FALSE)
  for(t in 1:M){
    H1[[t]] <- h1
    H2[[t]] <- h2
    X[[t]] <- muS_X + muT_X[t]*H1[[t]]*H2[[t]]/5 + .5*sample_field(chol_Sigma_exp,plt=FALSE)
    Y[[t]] <- X[[t]]*(3/2+H1[[t]]*H2[[t]]) + H1[[t]]^2 + abs(H2[[t]])*sample_field(chol_Sigma_exp,plt=FALSE)
  }
  x.out <- c(x.out, unlist(X))
  y.out <- c(y.out, unlist(Y))
  h1.out <- c(h1.out, unlist(H1))
  h2.out <- c(h2.out, unlist(H2))
}


########### 
## plotting
###########

if(1){
df1 <- expand.grid(s1=1:N1, s2=1:N2, t=1:M, b = 1:B)
df1$s <- rep(rep(1:(N), M*B))
df1$X <- x.out
df1$Y <- y.out
df1$H1 <- h1.out
df1$H2 <- h2.out
df1 <- df1[order(df1$b,df1$s1,df1$s2,df1$t),]
m <- 25
df1 <- subset(df1, t <= m & b == 1)
n1 <- n2 <- 25
df1 <- subset(df1, s1 <= n1 & s2 <= n2)

bhat.plt <- foreach(i = 1:(n1*n2), .combine = "rbind") %do% {
  print(paste(i, "out of ", n1*n2))
  xx <- df1$X[(m*(i-1)+1):((m*(i-1)+m))]
  yy <- df1$Y[(m*(i-1)+1):((m*(i-1)+m))]
  coefficients(lm(yy~xx))
}

b.causal.plt <- colMeans(bhat.plt)

linedata <- data.frame(intercept = bhat.plt[,1], 
                       slope = bhat.plt[,2])

p <- ggplot(df1, aes(X,Y)) + geom_point(alpha=.7, aes(col=H1*H2)) + 
  scale_color_gradient(name = expression(paste(bar(H)," * ",tilde(H))), low="black", high="gray") + 
  geom_abline(data=linedata, aes(slope=slope, intercept=intercept), col = "blue", alpha = .15, size = .5) + 
  geom_abline(intercept = mean(linedata$intercept), slope = mean(linedata$slope), col = "blue", size=1.5) + 
  theme(legend.position = c(.2,.7), 
        text = element_text(size=16)) + 
  geom_smooth(col = "red", se = FALSE, size = 1.5) + 
  geom_abline(intercept = b0[1], slope = b0[1], col = "green", size = 1.5, lty=2)
p

pdf("../figures/CSTM_example_XY.pdf", width=4*1.4, height = 4*1.3)
p
dev.off()

## X
pX.1 <- ggplot(subset(df1, t <=2), aes(s1, s2, fill=X,col=X)) + 
  geom_tile(size=.01) +
  scale_fill_distiller(palette = 'Spectral') +
  scale_color_distiller(palette = 'Spectral') +
  xlab("") + ylab("") +
  #coord_fixed() +
  facet_grid(t~.) +
  #ggtitle(expression(bold(X))) + 
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        strip.text = element_blank(),
        plot.title = element_text(hjust=.5), 
        plot.margin = unit(c(1,1,1,1), "mm"))
pX.1

pX.2 <- ggplot(subset(df1, t == m), aes(s1, s2, fill=X, col=X)) + 
  geom_tile(size=.01) +
  scale_fill_distiller(palette = 'Spectral') +
  scale_color_distiller(palette = 'Spectral') +
  xlab("") + ylab("") +
  #coord_fixed() +
  facet_grid(t~.) +
  #ggtitle(expression(bold(X))) + 
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        strip.text = element_blank(),
        plot.title = element_text(hjust=.5), 
        plot.margin = unit(c(1,1,1,1), "mm"))
pX.2

## Y 

pY.1 <- ggplot(subset(df1, t <=2), aes(s1, s2, fill=Y, col=Y)) + 
  geom_tile(size=.01) +
  scale_fill_distiller(palette = 'Spectral') +
  scale_color_distiller(palette = 'Spectral') +
  xlab("") + ylab("") +
  #coord_fixed() +
  facet_grid(t~.) +
  #ggtitle(expression(bold(X))) + 
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        strip.text = element_blank(),
        plot.title = element_text(hjust=.5), 
        plot.margin = unit(c(1,1,1,1), "mm"))
pY.1

pY.2 <- ggplot(subset(df1, t == m), aes(s1, s2, fill=Y, col=Y)) + 
  geom_tile(size=.01) +
  scale_fill_distiller(palette = 'Spectral') +
  scale_color_distiller(palette = 'Spectral') +
  xlab("") + ylab("") +
  #coord_fixed() +
  facet_grid(t~.) +
  #ggtitle(expression(bold(X))) + 
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        strip.text = element_blank(),
        plot.title = element_text(hjust=.5), 
        plot.margin = unit(c(1,1,1,1), "mm"))
pY.2


## H1
pH1.1 <- ggplot(subset(df1, t <=2), aes(s1, s2, fill=H1, col=H1)) + 
  geom_tile(size=.01) +
  scale_fill_distiller(palette = 'Spectral') +
  scale_color_distiller(palette = 'Spectral') +
  xlab("") + ylab("") +
  #coord_fixed() +
  facet_grid(t~.) +
  #ggtitle(expression(bold(X))) + 
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        strip.text = element_blank(),
        plot.title = element_text(hjust=.5), 
        plot.margin = unit(c(1,1,1,1), "mm"))
pH1.1

pH1.2 <- ggplot(subset(df1, t == m), aes(s1, s2, fill=H1, col=H1)) + 
  geom_tile(size=.01) +
  scale_fill_distiller(palette = 'Spectral') +
  scale_color_distiller(palette = 'Spectral') +
  xlab("") + ylab("") +
  #coord_fixed() +
  facet_grid(t~.) +
  #ggtitle(expression(bold(X))) + 
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        strip.text = element_blank(),
        plot.title = element_text(hjust=.5), 
        plot.margin = unit(c(1,1,1,1), "mm"))
pH1.2

## H2
pH2.1 <- ggplot(subset(df1, t <=2), aes(s1, s2, fill=H2, col=H2)) + 
  geom_tile(size=.01) +
  scale_fill_distiller(palette = 'Spectral') +
  scale_color_distiller(palette = 'Spectral') +
  xlab("") + ylab("") +
  #coord_fixed() +
  facet_grid(t~.) +
  #ggtitle(expression(bold(X))) + 
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        strip.text = element_blank(),
        plot.title = element_text(hjust=.5), 
        plot.margin = unit(c(1,1,1,1), "mm"))
pH2.1

pH2.2 <- ggplot(subset(df1, t == m), aes(s1, s2, fill=H2, col=H2)) + 
  geom_tile(size=.01) +
  scale_fill_distiller(palette = 'Spectral') +
  scale_color_distiller(palette = 'Spectral') +
  xlab("") + ylab("") +
  #coord_fixed() +
  facet_grid(t~.) +
  #ggtitle(expression(bold(X))) + 
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        strip.text = element_blank(),
        plot.title = element_text(hjust=.5), 
        plot.margin = unit(c(1,1,1,1), "mm"))
pH2.2



pdf("../figures/CSTM_maps_1.pdf", width=10, height = 4.4)
grid.arrange(pH1.1, pH2.1, pX.1, pY.1, nrow=1)
dev.off()

pdf("../figures/CSTM_maps_2.pdf", width=10, height = 2.3)
grid.arrange(pH1.2, pH2.2, pX.2, pY.2, nrow=1)
dev.off()

}

########################
## verifying convergence
########################

df <- expand.grid(s1=1:N1, s2=1:N2, t=1:M, b = 1:B)
df$s <- rep(rep(1:(N), M*B))
df$x <- x.out
df$y <- y.out
df <- df[order(df$b,df$s1,df$s2,df$t),]
x <- df$x
y <- df$y

# write.table(df, "CSTM_example_sim.txt", quote = FALSE)
# df <- read.table("CSTM_example_sim.txt", header = TRUE)

bhat <- foreach(i = 1:(B*N), .combine = "rbind") %do% {
  print(paste(i, "out of ", B*N))
  foreach(j = 1:(M/10), .combine = "rbind") %do% {
    xx <- x[(M*(i-1)+1):((M*(i-1)+j*10))]
    yy <- y[(M*(i-1)+1):((M*(i-1)+j*10))]
    coefficients(lm(yy~xx))
  }
}


df.sub <- subset(df, t %in% seq(10,M,10))
df.sub$b0hat <- bhat[,1]
df.sub$b1hat <- bhat[,2]
df.sub <- df.sub[order(df.sub$b,df.sub$t,df.sub$s1,df.sub$s2),]

write.table(df.sub, "CSTM_example_sim.txt", quote = FALSE)
df.sub <- read.table("CSTM_example_sim.txt", header = TRUE)


b0hat <- df.sub$b0hat
b1hat <- df.sub$b1hat

b.causal.hat <- foreach(i = 1:(B*M/10), .combine = "rbind") %do% {
  print(paste(i, "out of ", B*M/10))
  foreach(j = 1:(N1/10), .combine = "rbind") %do% {
    c(mean(b0hat[(N*(i-1)+1):(N*(i-1)+N2*j*10)],na.rm=1), 
      mean(b1hat[(N*(i-1)+1):(N*(i-1)+N2*j*10)],na.rm=1))
  }
}

df.avgce <- subset(df.sub, s2 == 1 & s1 %in% seq(10,N1,10))
df.avgce$b.causal0.hat <- b.causal.hat[,1]
df.avgce$b.causal1.hat <- b.causal.hat[,2]
df.avgce$err <- sqrt((df.avgce$b.causal0.hat-b0)^2 + (df.avgce$b.causal1.hat-b1)^2)

write.table(df.avgce, "CSTM_example_bhat_frame.txt", quote = FALSE)
df.avgce <- read.table("CSTM_example_bhat_frame.txt", header = TRUE)

eps <- .2
df.quant <- aggregate(err ~ s1 + t, df.avgce, function(e) mean(e > eps))
df.quant$n <- df.quant$s1*N2

pp <- ggplot(df.quant, aes(n,t,fill=err,col=err)) + 
  geom_tile(size=.01) +
  # scale_fill_distiller(name = "err. prob.", palette = 'Spectral', type = "seq") +
  scale_fill_distiller(name = "error prob.", palette = "YlOrRd", direction = 1) + 
  scale_color_distiller(name = "error prob.", palette = "YlOrRd", direction = 1) + 
  # scale_fill_gradient2(name = "err. prob.", low = "green", mid = "yellow", high = "red", midpoint = .35) +
  theme(text = element_text(size=16)) + 
  xlab("number of spatial locations") + 
  ylab("number of temporal instances") 
pp

pdf("../figures/CSTM_example_convergence.pdf", width=4*1.7, height = 4*1.3)
pp
dev.off()
