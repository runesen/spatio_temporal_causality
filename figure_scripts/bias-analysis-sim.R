setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(foreach)
library(gridExtra)
library(reshape2)
library(EnvStats)
theme_set(theme_bw())


dat <- read.table("data_xy_colombia_20200616.txt", header = TRUE)
dat <- dat[order(dat$PolygonID),]
n <- nrow(dat)
nt <-  length(unique(dat$Year))
ns <- length(unique(dat$PolygonID))

X <- dat$nr_fatalities>0 # binary conflict indicator
Y <- dat$FL_km # forestloss
H <- dat$RoadDist


########### method statistics ###########

# to compute the lscm method-statistic, we can disregard many observations
Xmat <- matrix(X, nrow=nt, ncol=ns)
Xmatt <- matrix(rep(colSums(Xmat),each=nt), nrow=nt, ncol=ns)
include <- c(Xmatt)>0 & c(Xmatt) < nt # true for all locations containing conflicts as well as no conflicts


## method-statistic using difference in mean
ts.lscm <- function(X, Y){
  
  X <- X[include]
  Y <- Y[include]
  ns <- length(X)/nt
  
  dfa <- foreach(s=1:ns, .combine = "c") %do% {
    w <- (s-1)*nt + 1:nt
    Xs <- X[w]
    Ys <- Y[w]
    w0 <- which(!Xs & !is.na(Ys))
    w1 <- which(Xs & !is.na(Ys))
    mean(Ys[w1])-mean(Ys[w0])
  }
  mean(dfa,na.rm=1)
} 

## only confounder: H
n.H.quant <- 100
qseq <- quantile(H, prob = seq(0,1-1/n.H.quant,1/n.H.quant)) # quantiles
w <- sapply(1:n, function(i) max(which(H[i]-qseq >=0))) # assign each data point to a quantile index
H.quantile.groups <- lapply(1:length(qseq), function(i) which(w==i)) # lists with indices of data points belonging to respective quantile group
n.confl.H.quantile.groups <- sapply(H.quantile.groups, function(i) sum(X[i])) # removing quantile groups without conflict
H.quantile.groups <- H.quantile.groups[n.confl.H.quantile.groups>0]
n.H.quant <- length(H.quantile.groups)

ts.H <- function(X,Y){
  # difference in sample averages per quantile group
  t <- sapply(1:n.H.quant, function(i){
    ind.qi <- H.quantile.groups[[i]]
    confl.qi <- intersect(which(X), ind.qi)
    noconfl.qi <- intersect(which(!X), ind.qi)
    mean(Y[confl.qi],na.rm=1) - mean(Y[noconfl.qi],na.rm=1)
  })
  # remove possible NAs
  w.not.NA <- which(!is.na(t))
  t <- t[w.not.NA]
  H.quantile.groups <- H.quantile.groups[w.not.NA]
  weights <- sapply(H.quantile.groups, length) # empirical distribution of discretized population density
  weights <- weights / sum(weights)
  sum(t*weights)
}

## no confounder
ts.noconf <- function(X,Y){
  mean(Y[X],na.rm=1)-mean(Y[!X],na.rm=1)
}



################### 
## results on actual data
################### 
(ts.noconf.data <- ts.noconf(X,Y)) # 0.073
(ts.H.data <- ts.H(X,Y)) # 0.039
(ts.lscm.data <- ts.lscm(X,Y)) # -0.0180
####################### power analysis #######################



# function that simulates pseudo-data sets where the true method statistic (the difference in avg. causal effects between X=1 and X=0) is b
Ysim.fun <- function(b, a){
  ind.res <- c(sapply(1:ns, function(s) sample(nt*(s-1)+1:nt, replace = TRUE)))
  b*X + a*log(H) + Y[ind.res] # simulated response
}

hist(log(H))
hist(log(H)[include])


################### 
## bias analysis
################### 
n.b <- 5
(a.seq <- seq(-.02, .02, length.out = n.b))
(b.seq <- seq(-.2, .2, length.out = n.b))
n.sim <- 100
method.vec <- c("conf", "adjustH", "lscm")
ts.lst <- list(ts.noconf, ts.H, ts.lscm)
n.method <- length(method.vec)
s <- 0
out <- NULL
for(i in 1:n.sim){
  ts.true.out <- ts.hat.out <- NULL
  for(k in 1:n.method){
    method <- method.vec[k]
    ts <- ts.lst[[k]]
    for(b in b.seq){
      for(a in a.seq){
        s <- s+1
        set.seed(s)
        print(paste("method = ", method, ", SIM = ", i, ", b = ", b, ", a = ", a))
        Ysim <- Ysim.fun(b,a)
        ts.hat <- ts(X, Ysim)
        ts.true <- b #+ a*mean(log(H))
        # print(paste("ts.true = ", ts.true, ", ts.hat = ", ts.hat))
        ts.true.out <- c(ts.true.out, ts.true)
        ts.hat.out <- c(ts.hat.out, ts.hat)
      }
    }
  }
  tmp <- data.frame(method = rep(method.vec, each = n.b*n.b), 
                    b = rep(rep(b.seq, n.method), each=n.b),
                    a = rep(a.seq, n.b*n.method),
                    ts.hat = ts.hat.out, 
                    ts.true = ts.true.out, 
                    bias = ts.hat.out - ts.true.out,
                    sim = i)
  out <- rbind(out, tmp)
  write.table(out, "2021-02-19_bias-analisys.txt", quote = FALSE)
}