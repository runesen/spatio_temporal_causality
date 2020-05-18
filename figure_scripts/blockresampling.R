###################################
# block resampling procedure
##################################"
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(foreach)
library(gridExtra)
library(reshape2)
theme_set(theme_bw())


# dat1 <- read.table("data_xy_colombia_20191219.txt", header = TRUE)
# dat1 <- read.table("data_xy_colombia_20200316.txt", header = TRUE)
dat1 <- read.table("data_xy_colombia_20200327.txt", header = TRUE)
w <- which(dat1$xFor2000_ge25==0)
dat1 <- dat1[-w,]
n <- nrow(dat1)
nt <- ny <- length(unique(dat1$Year))
ns <- length(unique(dat1$PolygonID))

c <- dat1$nr_fatalities>0 # binary conflict indicator
f <- dat1$xFor_ge25 # forest
fl <- dat1$FL_km # forestloss


Ts <- function(c,fl){ # here fl must be lagged, i.e., fl[1] is f[2]-f[1]
  if(var(c)==0) return(NA)
  else{
    w0 <- which(!c & !is.na(fl))
    w1 <- which(c & !is.na(fl))
    mean(fl[w1])-mean(fl[w0])
  }
}

cmat <- matrix(c, nrow=nt, ncol=ns)
cmatt <- matrix(rep(colSums(cmat),each=nt), nrow=nt, ncol=ns)
ind <- c(cmatt)>=1 & c(cmatt)<nt # true for all datapoints from conflict locations


Tpixel <- function(c,fl){ 
  
  c <- c[ind]
  fl <- fl[ind]
  ns <- length(c)/nt
  
  dfa <- foreach(s=1:ns, .combine = "c") %do% {
    w <- (s-1)*nt + 1:nt
    cs <- c[w]
    fls <- fl[w]
    # Ts(cs[-nt],fls[-1])
    Ts(cs,fls)
  }
  dfa
}


tpixelvec <- Tpixel(c,fl)
(tpixel <- mean(tpixelvec ,na.rm=1))

nblocks <- 6
ind.lst <- list(1:3,4:6,7:9,10:12,13:15,16:19)

B <- 999
Tpixelboot <- foreach(b=0:B, .combine="rbind") %do% {
  
  print(paste(b, " out of ", B))
  set.seed(b)
  if(b == 0){ # actual data 
    zz <- 1:(nt*ns)
  } else{ # resampled data
    z <- unlist(ind.lst[sample(1:nblocks, size = nblocks, replace=FALSE)]) # blockresampling
    # shft <- sample(0:2, size=1)
    shft <- 0
    z <- c(tail(z,shft), head(z,nt-shft)) # random shift
    zz <- nt*rep(0:(ns-1), each=nt) + rep(z, ns)
  }
  
  flb <- fl[zz]
  
  mean(Tpixel(c,flb),na.rm=1)
}

hist(Tpixelboot[-1], breaks = 20)
abline(v=Tpixelboot[1], col="red")

bootframe <- data.frame(t = Tpixelboot,
                        b = 0:B, 
                        type = "blockresampling")
rownames(bootframe) <- NULL


# write.table(bootframe, "blockresampling_data_new.txt", quote = FALSE)
# bootframe <- read.table("blockresampling_data_new.txt", header = TRUE)


## p-value
t <- bootframe$t
pv <- min(1,2*min((1+sum(t[-1] <= t[1]))/(1+B), (1+sum(t[-1] >= t[1]))/(1+B)))
t[1] # -0.0180
pv # 0.892
