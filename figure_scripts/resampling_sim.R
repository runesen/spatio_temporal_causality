###################################
# calculating test statistics for different
# null hypothesis alongside with permutations
##################################"
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(foreach)
library(gridExtra)
library(reshape2)
theme_set(theme_bw())


dat <- read.table("data_xy_colombia_20200327.txt", header = TRUE)
dat <- dat[order(dat$PolygonID),]
n <- nrow(dat)
nt <-  length(unique(dat$Year))
ns <- length(unique(dat$PolygonID))

c <- dat$nr_fatalities>0 # binary conflict indicator
f <- dat$xFor_ge25 # forest
fl <- dat$FL_km # forestloss
popdens <- dat$PopDens
roaddist <- dat$RoadDist


############ absolute forest loss ############ 

# average deforestation rate
mean(fl,na.rm=1) # 0.193
# average deforestation at conflict
mean(fl[c],na.rm=1) # 0.266
# average deforestation NOT at conflict
mean(fl[!c],na.rm=1) # 0.193
# t-test
t.test(fl~c)$p.value # 8.416e-10



########### test statistics ###########

##########
## no confounder
##########
ts.noconf <- function(c,fl){
  mean(fl[c],na.rm=1)-mean(fl[!c],na.rm=1)
}

##########
## only confounder: population density
##########
n.popdens.quant <- 100
qseq <- quantile(popdens, prob = seq(0,1-1/n.popdens.quant,1/n.popdens.quant)) # quantiles
w <- sapply(1:n, function(i) max(which(popdens[i]-qseq >=0))) # assign each data point to a quantile index
popdens.quantile.groups <- lapply(1:n.popdens.quant, function(i) which(w==i)) # lists with indices of data points belonging to respective quantile group
n.confl.popdens.quantile.groups <- sapply(popdens.quantile.groups, function(i) sum(c[i])) # removing quantile groups without conflict
popdens.quantile.groups <- popdens.quantile.groups[n.confl.popdens.quantile.groups>0]
n.popdens.quant <- length(popdens.quantile.groups)
zz <- unlist(popdens.quantile.groups)

ts.popdens <- function(c,fl){
  # difference in sample averages per quantile group
  t <- sapply(1:n.popdens.quant, function(i){
    ind.qi <- popdens.quantile.groups[[i]]
    confl.qi <- intersect(which(c), ind.qi)
    noconfl.qi <- intersect(which(!c), ind.qi)
    mean(fl[confl.qi],na.rm=1) - mean(fl[noconfl.qi],na.rm=1)
  })
  # remove possible NAs
  w.not.NA <- which(!is.na(t))
  t <- t[w.not.NA]
  popdens.quantile.groups <- popdens.quantile.groups[w.not.NA]
  weights <- sapply(popdens.quantile.groups, length) # empirical distribution of discretized population density
  weights <- weights / sum(weights)
  sum(t*weights)
}

##########
## only confounder: distance to road
##########
n.roaddist.quant <- 100
qseq <- quantile(roaddist, prob = seq(0,1-1/n.roaddist.quant,1/n.roaddist.quant)) # quantiles
w <- sapply(1:n, function(i) max(which(roaddist[i]-qseq >=0))) # assign each data point to a quantile index
roaddist.quantile.groups <- lapply(1:length(qseq), function(i) which(w==i)) # lists with indices of data points belonging to respective quantile group
n.confl.roaddist.quantile.groups <- sapply(roaddist.quantile.groups, function(i) sum(c[i])) # removing quantile groups without conflict
roaddist.quantile.groups <- roaddist.quantile.groups[n.confl.roaddist.quantile.groups>0]
n.roaddist.quant <- length(roaddist.quantile.groups)

ts.roaddist <- function(c,fl){
  # difference in sample averages per quantile group
  t <- sapply(1:n.roaddist.quant, function(i){
    ind.qi <- roaddist.quantile.groups[[i]]
    confl.qi <- intersect(which(c), ind.qi)
    noconfl.qi <- intersect(which(!c), ind.qi)
    mean(fl[confl.qi],na.rm=1) - mean(fl[noconfl.qi],na.rm=1)
  })
  # remove possible NAs
  w.not.NA <- which(!is.na(t))
  t <- t[w.not.NA]
  roaddist.quantile.groups <- roaddist.quantile.groups[w.not.NA]
  weights <- sapply(roaddist.quantile.groups, length) # empirical distribution of discretized population density
  weights <- weights / sum(weights)
  sum(t*weights)
}

##########
## LSCM
##########
cmat <- matrix(c, nrow=nt, ncol=ns)
cmatt <- matrix(rep(colSums(cmat),each=nt), nrow=nt, ncol=ns)
ind <- c(cmatt)>0 & c(cmatt) < nt # true for all locations containing conflicts as well as no conflicts

ts.lscm <- function(c, fl){
  # removing irrelevant data 
  c <- c[ind] 
  fl <- fl[ind]
  ns.ind <- length(c)/nt
  
  dfa <- foreach(s=1:ns.ind, .combine = "c") %do% {
    w <- (s-1)*nt + 1:nt
    cs <- c[w]
    fls <- fl[w]
    w0 <- which(!cs & !is.na(fls))
    w1 <- which(cs & !is.na(fls))
    mean(fls[w1])-mean(fls[w0])
  }
  mean(dfa,na.rm=1)
}

########################
(ts.noconf.data <- ts.noconf(c,fl)) # 0.073
(ts.popdens.data <- ts.popdens(c,fl)) # 0.038
(ts.roaddist.data <- ts.roaddist(c,fl)) # 0.039
(ts.lscm.data <- ts.lscm(c,fl)) # -0.0180
########################

############ resampling procedures ############

##########
## no confounder
##########
res.noconf <- function(fl){
  # random resampling
  ind.res <- sample(1:n, replace = FALSE)
  fl[ind.res]
}

##########
## only conf: population density
##########
res.popdens <- function(fl){
  # resampling within quantile groups of popdens
  fl.res <- fl
  for(i in 1:n.popdens.quant){
    ind.qi <- popdens.quantile.groups[[i]]
    ind.qi.res <- sample(ind.qi, replace = FALSE)
    fl.res[ind.qi] <- fl[ind.qi.res]
  }
  fl.res
}


##########
## only conf: distance to road
##########
res.roaddist <- function(fl){
  # resampling within quantile groups of roaddist
  for(i in 1:n.roaddist.quant){
    ind.qi <- roaddist.quantile.groups[[i]]
    ind.qi.res <- sample(ind.qi, replace = FALSE)
    fl[ind.qi] <- fl[ind.qi.res]
  }
  fl
}

##########
## LSCM
##########
res.lscm <- function(fl){
  # resampling within the same location (i.e., values of the hidden confounders)
  time.res <- sample(1:nt, replace=FALSE)
  ind.res <- nt*rep(0:(ns-1), each=nt) + rep(time.res, ns)
  fl[ind.res]
}

########################## resampling test statistics ##########################
B <- 999
ts.resampler <- function(c, fl, test.stat, resampler){
  foreach(b = 1:B, .combine = "c") %do% {
    print(paste(b, " out of ", B))
    set.seed(b)
    fl.res <- resampler(fl)
    test.stat(c,fl.res)
  }
}

ts.noconf.res <- ts.resampler(c,fl,ts.noconf, res.noconf)
ts.popdens.res <- ts.resampler(c,fl,ts.popdens, res.popdens)
ts.roaddist.res <- ts.resampler(c,fl,ts.roaddist, res.roaddist)
ts.lscm.res <- ts.resampler(c,fl,ts.lscm, res.lscm)

xmin <- min(ts.noconf.data, 
            ts.popdens.data, 
            ts.roaddist.data, 
            ts.lscm.data, 
            ts.noconf.res, 
            ts.popdens.res, 
            ts.roaddist.res, 
            ts.lscm.res)
xmax <- max(ts.noconf.data, 
            ts.popdens.data, 
            ts.roaddist.data, 
            ts.lscm.data, 
            ts.noconf.res, 
            ts.popdens.res, 
            ts.roaddist.res, 
            ts.lscm.res)

par(mfrow=c(1,4))
hist(ts.noconf.res, breaks = 20, xlim=c(xmin, xmax))
abline(v=ts.noconf.data, col="red",lwd=3)
hist(ts.popdens.res, breaks = 20, xlim=c(xmin, xmax))
abline(v=ts.popdens.data, col="red",lwd=3)
hist(ts.roaddist.res, breaks = 20, xlim=c(xmin, xmax))
abline(v=ts.roaddist.data, col="red",lwd=3)
hist(ts.lscm.res, breaks = 20, xlim=c(xmin, xmax))
abline(v=ts.lscm.data, col="red",lwd=3)

res.frame <- data.frame(ts = c(ts.noconf.data, ts.noconf.res, 
                              ts.popdens.data, ts.popdens.res, 
                              ts.roaddist.data, ts.roaddist.res, 
                              ts.lscm.data, ts.lscm.res), 
                        b = rep(0:B, 4), 
                        model = rep(c("noconf", "popdens", "roaddist", "lscm"), each = B+1)
                        )

write.table(res.frame, "resampling_data.txt", quote = FALSE)
