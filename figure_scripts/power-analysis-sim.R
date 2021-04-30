setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(foreach)
library(gridExtra)
library(reshape2)
library(EnvStats)
source("resampling-test.R")
theme_set(theme_bw())

dat <- read.table("data_xy_colombia_20200616.txt", header = TRUE)
dat <- dat[order(dat$PolygonID),]
n <- nrow(dat)
nt <-  length(unique(dat$Year))
ns <- length(unique(dat$PolygonID))

X <- dat$nr_fatalities>0 # binary conflict indicator
Y <- dat$FL_km # forestloss
H <- dat$RoadDist


########### test statistics ###########

# to compute the lscm test-statistic, we can disregard many observations
Xmat <- matrix(X, nrow=nt, ncol=ns)
Xmatt <- matrix(rep(colSums(Xmat),each=nt), nrow=nt, ncol=ns)
include <- c(Xmatt)>0 & c(Xmatt) < nt # true for all locations containing conflicts as well as no conflicts


## test-statistic using difference in mean
ts.lscm <- function(X, Y=Ysim.res){
  
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
    # mean(Ys[w1]>0)-mean(Ys[w0]>0)
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


############ resampling procedures ############

##########
## LSCM: resampling along time axis
##########
res.procedure.lscm <- function(Y){
  # resampling within the same location (i.e., values of the hidden confounders)
  time.res <- sample(1:nt, replace=FALSE)
  ind.res <- nt*rep(0:(ns-1), each=nt) + rep(time.res, ns)
  Y[ind.res]
}

Ysim.res <- res.procedure.lscm(Ysim)


##########
## no confounder
##########
res.procedure.noconf <- function(Y){
  # random resampling
  ind.res <- sample(1:n, replace = FALSE)
  Y[ind.res]
}

##########
## only conf: distance to road
##########
res.procedure.H <- function(Y){
  # resampling within quantile groups of roaddist
  for(i in 1:n.H.quant){
    ind.qi <- H.quantile.groups[[i]]
    ind.qi.res <- sample(ind.qi, replace = FALSE)
    Y[ind.qi] <- Y[ind.qi.res]
  }
  Y
}


####################### resampling test #######################
## computing resampled test statistics
ts.resampler <- function(X, Y, test.stat, res.procedure, B=99){
  foreach(b = 1:B, .combine = "c") %do% {
    if(b %% 10 == 0) print(paste(b, " out of ", B))
    #set.seed(b)
    Y.res <- res.procedure(Y)
    test.stat(X,Y.res)
  }
}

## p-value for two-sided test
pvfun <- function(ts.data, ts.res, na.rm){
  ts.res <- ts.res[!is.na(ts.res)]
  B <- length(ts.res)
  min(1,2*min((1+sum(ts.res <= ts.data))/(1+B), (1+sum(ts.res >= ts.data))/(1+B)))
} 



################### 
## results on actual data
################### 
(ts.noconf.data <- ts.noconf(X,Y)) # 0.073
(ts.H.data <- ts.H(X,Y)) # 0.039
(ts.lscm.data <- ts.lscm(X,Y)) # -0.0180

ts.noconf.res <- ts.resampler(X,Y,ts.noconf, res.procedure.noconf)
ts.H.res <- ts.resampler(X,Y,ts.H, res.procedure.H)
ts.lscm.res <- ts.resampler(X,Y,ts.lscm, res.procedure.lscm)


(pv.noconf <- pvfun(ts.noconf.data, ts.noconf.res, na.rm=TRUE))
(pv.H <- pvfun(ts.H.data, ts.H.res, na.rm=TRUE))
(pv.lscm <- pvfun(ts.lscm.data, ts.lscm.res, na.rm=TRUE))




####################### power analysis #######################

# procedure for generating responses
Ysim.fun <- function(b){
  ind.res <- c(sapply(1:ns, function(s) sample(nt*(s-1)+1:nt, replace = TRUE)))
  b*X + Y[ind.res] # simulated response
}


################### 
## power analysis
################### 
B <- 99
n.b <- 5
(b.seq <- seq(0, .08, length.out = n.b))
n.sim <- 100
test.vec <- c("t-test", "resampling-conconf", "resampling-adjustH", "resampling-lscm")
ts.lst <- list(ts.noconf, ts.noconf, ts.H, ts.lscm)
res.procedure.lst <- list(NA, res.procedure.noconf, res.procedure.H, res.procedure.lscm)
n.ts <- length(test.vec)
s <- 0
out <- NULL
for(i in 1:n.sim){
  pv.out <- NULL
  ts.hat.out <- NULL
  for(k in 1:n.ts){
    test <- test.vec[k]
    res.procedure <- res.procedure.lst[[k]]
    ts <- ts.lst[[k]]
    for(b in b.seq){
      s <- s+1
      set.seed(s)
      print(paste("test = ", test, ", SIM = ", i, 
                  
                  ", b = ", b))
      Ysim <- Ysim.fun(b)
      ts.hat <- ts(X, Ysim)
      if(test == "t-test"){
        pv <- t.test(Ysim ~ X)$p.value
      } else{
        ts.res <- ts.resampler(X, Ysim, ts, res.procedure, B)
        pv <- pvfun(ts.hat, ts.res)
      }
      print(paste("pval = ", pv))
      ts.hat.out <- c(ts.hat.out, ts.hat)
      pv.out <- c(pv.out, pv)
    }
  }
  tmp <- data.frame(test = rep(test.vec, each = n.b), 
                    b = rep(b.seq, n.ts), 
                    pv = pv.out, 
                    ts.hat = ts.hat.out, 
                    sim = i)
  out <- rbind(out, tmp)
  write.table(out, "2021-02-15_power-analisys.txt", quote = FALSE)
}
