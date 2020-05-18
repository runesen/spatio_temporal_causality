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


df <- read.table("data_xy_colombia_20200327.txt", header = TRUE)
n <- nrow(df)
nt <- length(unique(df$Year))
ns <- length(unique(df$PolygonID))
d.pixel <- 10000
fact <- 10 # sidelength of quadratic blocks that will be permuted together
lat.rd <- min(df$lat) + d.pixel*fact*floor((df$lat-min(df$lat)) / (fact*d.pixel)+.0001) + d.pixel/2 # rounded lat
lon.rd <- min(df$lon) + d.pixel*fact*floor((df$lon-min(df$lon)) / (fact*d.pixel)+.0001) + d.pixel/2 # and lon
df$square <- interaction(lat.rd, lon.rd) # square indicator
df$inner <- table(df$square)[df$square]==nt*fact^2 # indicates whether square is full of data points


# tmp <- data.frame(lat = rep(unique(lat.rd), each=length(unique(lon.rd))), 
#                   lon = rep(unique(lon.rd), length(unique(lat.rd))))
# tmp <- tmp[order(tmp$lon, tmp$lat),]
# tmp$count <- as.numeric(table(df$square))

#####################################
## visualization of spatial blocks
#####################################

df.sq <- aggregate(lat~square, subset(df,inner), mean)
df.sq$lon <- aggregate(lon~square, subset(df,inner), mean)$lon
df.sq$square <- 1:nrow(df.sq)

write.table(df.sq, "spatial-blocks.txt", quote=FALSE)

p <- ggplot() + 
  geom_rect(data=subset(df,Year==2000), aes(xmin=lon-d.pixel/2, xmax=lon+d.pixel/2, ymin=lat-d.pixel/2, ymax=lat+d.pixel/2), alpha=.5, size=.05) + 
  geom_rect(data=df.sq, aes(xmin=lon-fact*d.pixel/2, xmax=lon+fact*d.pixel/2, ymin=lat-fact*d.pixel/2, ymax=lat+fact*d.pixel/2), alpha=.5, col="blue") + 
  # geom_tile(data=subset(df,Year==2000), aes(lon, lat), alpha=.5) + 
  # geom_tile(data=df.sq, aes(lon,lat), alpha=1, col = "black") + 
  scale_alpha_manual(limits = c(TRUE, FALSE), values = c(1,.5)) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none", 
        plot.margin = margin(t=2, r=15, b=2, l=2, unit="mm")) + 
  ggtitle("block structure for spatial permutations") + 
  xlab("") + ylab("")
p

#####################################
## permutation tests
#####################################

c <- df$nr_fatalities>0
fl <- df$FL_km
sq <- df$square
inner <- df$inner

# data ordered such that all inner squares come first
c <- c[order(!inner,sq)] 
fl <- fl[order(!inner,sq)]

num.blocks <- sum(inner)/fact^2
w.inner <- 1:(num.blocks*fact^2)
w.outer <- (num.blocks*fact^2+1):n

# test statistic
ts <- function(c,fl){
  mean(fl[c],na.rm=1)-mean(fl[!c],na.rm=1)
}

ts.data <- ts(c,fl)

B <- 999
#####################################
# spatial blockresampling of FL
# also permuting boundary points
#####################################


res.random <- function(fl){
  fl.res <- fl
  rand.ind <- sample(1:num.blocks, size = num.blocks, replace = FALSE)
  res.ind.inner <- c(sapply(rand.ind, function(i) (i-1)*fact^2+1:(fact^2)))
  res.ind.outer <- sample(w.outer, size = length(w.outer), replace=FALSE)
  fl.res[w.inner] <- fl[res.ind.inner]
  fl.res[w.outer] <- fl[res.ind.outer]
  fl.res
}

ts.res.block <- foreach(b=1:B, .combine = "c") %do% {
  print(paste("b = ",b))
  set.seed(b)
  fl.res <- res.random(fl)
  ts(c,fl.res)
}
# loading result from the resampling test performed in the script resampling_sim.R
ts.res.random <- subset(read.table("resampling_data_new.txt", header = TRUE), test=="marginal" & b > 0)$t

df.plt <- data.frame(ts = c(ts.res.block, ts.res.random), 
                     method = rep(c("block permutation", "random permutation"), each=B))

p.res <- ggplot(df.plt, aes(x=ts)) + geom_histogram(color="gray", alpha=.7) + 
  facet_grid(method~.) + 
  geom_vline(xintercept = ts.data, col="red", size = 1) + 
  xlab("test statistic") + 
  ggtitle("resampling tests") + 
  theme(#text = element_text(size=15), 
        plot.title = element_text(hjust=.5), 
        plot.margin = margin(t=2, r=2, b=2, l=15, unit="mm"))

pdf("../figures/spatial_block_resampling_all.pdf", width = 10, height = 6)
grid.arrange(p, p.res,ncol=2)
dev.off()

(pv.block <- min(1,2*min((1+sum(ts.res.block <= ts.data))/(1+B), (1+sum(ts.res.block >= ts.data))/(1+B)))) # 0.008
(pv.random <- min(1,2*min((1+sum(ts.res.random <= ts.data))/(1+B), (1+sum(ts.res.random >= ts.data))/(1+B)))) # 0.002

