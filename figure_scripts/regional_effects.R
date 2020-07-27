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

df <- read.table("data_xy_colombia_20200616.txt", header = T)
df.prov <- read.table("uniqueID_country-province.txt", header = TRUE, sep = "\t")
df.prov <- subset(df.prov, GID_0 == "COL")
ns <- length(unique(df$PolygonID))
nt <- length(unique(df$Year))
df <- df[order(df$PolygonID),]
df.prov <- df.prov[order(df.prov$UniqueID),]
df$province <- factor(rep(as.character(df.prov$NAME_1), each = nt))
# levels(df$province)[4]  <- "Atlantico"
# levels(df$province)[6]  <- "Boyaca"
# levels(df$province)[8]  <- "Caqueta"
# levels(df$province)[21]  <- "Narino"

df.agg <- aggregate(nr_fatalities ~ province, df, function(x) sum(x>0))
df.agg$lat <- aggregate(lat ~ province, df, mean)$lat
df.agg$lon <- aggregate(lon ~ province, df, mean)$lon
df.agg$risk <- df.agg$nr_fatalities / aggregate(PolygonID ~ province, df, length)$PolygonID

ggplot() + geom_tile(data=subset(df, Year == 2000), aes(lon, lat, fill=province)) + 
  geom_text(data=df.agg, aes(x=lon, y=lat, label=nr_fatalities))

w.c <- which(aggregate(nr_fatalities ~ PolygonID, df, function(x) any(x>0))$nr_fatalities)
poly.c <- unique(df$PolygonID)[w.c]
df.sub <- subset(df, PolygonID %in% poly.c)
ns.c <- length(poly.c)
c <- df.sub$nr_fatalities>0
fl <- df.sub$FL_km
prov <- as.character(subset(df.sub, Year == 2000)$province)

# computes diff in sample avg for a given pixel
Ts <- function(c=cs,fl=fls){ 
  if(var(c)==0) return(NA)
  else{
    w0 <- which(!c & !is.na(fl))
    w1 <- which(c & !is.na(fl))
    mean(fl[w1])-mean(fl[w0])
  }
}

test.stats <- function(c,fl){
  ts <- foreach(s=1:ns.c, .combine = "c") %do% {
    w <- (s-1)*nt + 1:nt
    cs <- c[w]
    fls <- fl[w]
    Ts(cs,fls)
  }
  aggregate(ts~prov, FUN=mean, na.rm=1)$ts
}

test.stats(c,fl)

n.p <- length(unique(prov))

B <- 999
test.stats.boot <- foreach(b=0:B, .combine="rbind") %do% {
  
  print(paste(b, " out of ", B))
  set.seed(b)
  if(b == 0){ # actual data 
    zz <- 1:(nt*ns.c)
  } else{ # resampled data
    z <- sample(1:nt, size = nt, replace=FALSE)
    zz <- nt*rep(0:(ns.c-1), each=nt) + rep(z, ns.c)
  }
  flb <- fl[zz]
  tt <- test.stats(c,flb)
  if(length(tt) != n.p) tt <- rep(NA, n.p)
  tt
}

df.test.stat <- as.data.frame(test.stats.boot)
colnames(df.test.stat) <- sort(unique(prov))
df.test.stat$b <- 0:B
df.test.stat.melt <- melt(df.test.stat, id.vars = "b")
colnames(df.test.stat.melt) <- c("b", "province", "t")

pvfun <- function(t, na.rm){
  t.data <- t[1]
  t.boot <- t[-1]
  B.non.NA <- sum(!is.na(t.boot))
  min(1,2*min((1+sum(t.boot <= t.data))/(1+B.non.NA), (1+sum(t.boot >= t.data))/(1+B.non.NA)))
} 

pvframe <- aggregate(t ~ province, df.test.stat.melt, FUN = pvfun)
colnames(pvframe)[2] <- "pval"
pvframe$t <- subset(df.test.stat.melt, b==0)$t

prov.NA <- setdiff(df.plt$province, pvframe$province)
pvframe.NA <- data.frame(province = prov.NA, pval = NA, t = NA)
pvframe <- rbind(pvframe, pvframe.NA)

write.table(pvframe, "regional_effects.txt", quote = FALSE, sep = "\t")
pvframe <- read.table("regional_effects.txt", header = TRUE, sep = "\t")

ggplot(subset(df.test.stat.melt, b > 0), aes(x=t)) + geom_histogram(color="gray", alpha=0.7) + 
  facet_wrap(.~province) + 
  geom_vline(data=pvframe, aes(xintercept=t, col=t), size = 1) + 
  scale_color_distiller(palette="RdYlBu") + 
  geom_text(data=pvframe, aes(x=0, y=B/2,label=round(pval,3)), size = 5) + 
  xlab("test statistic")


df.plt <- subset(df, Year == 2000)
df.plt$province <- as.character(df.plt$province)
df.plt <- df.plt[order(df.plt$province),]
pvframe <- pvframe[order(as.character(pvframe$province)),]
n.prov <- table(df.plt$province)
df.plt$t <- unlist(sapply(1:length(n.prov), function(i) rep(pvframe$t[i], n.prov[i])))
df.plt$pval <- unlist(sapply(1:length(n.prov), function(i) rep(pvframe$pval[i], n.prov[i])))
df.plt$effect <- factor(sign(df.plt$t), levels = c(-1,1), labels = c("positive", "negative"))
df.plt$significant <- factor(df.plt$pval <= 0.05, levels = c(TRUE, FALSE), labels = c("true", "false"))

df.agg1 <- aggregate(nr_fatalities ~ lat + lon, df, function(x) any(x>0))
df.agg1 <- subset(df.agg1, nr_fatalities)[,1:2]
write.table(df.agg1, "colombia_conflict_locations_2000-2018.txt", quote = FALSE)


pdf("../figures/regional_effects.pdf", width = 5.07, height = 5.8)
ggplot() + geom_tile(data=df.plt, aes(lon, lat, fill=effect, alpha=significant)) + 
  geom_point(data = df.agg1, aes(lon, lat), col = "red", shape = 4) + 
  geom_text(data=df.agg, aes(x=lon, y=lat, label=nr_fatalities)) +
  scale_fill_manual(values = c("#D55E00", "#009E73"), na.value = "gray") + 
  scale_alpha_manual(values = c(1,.5), na.value = .5) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_blank(), 
        axis.ticks = element_blank()) + 
  xlab("") + ylab("")
dev.off()


ggplot() + geom_tile(data=df.plt, aes(lon, lat, fill=t)) + 
  # geom_point(data = df.agg1, aes(lon, lat), col = "red", shape = 4) + 
  # geom_text(data=df.agg, aes(x=lon, y=lat, label=nr_fatalities)) +
  # scale_fill_manual(values = c("#D55E00", "#009E73"), na.value = "gray") + 
  scale_fill_distiller(palette = "RdYlGn") + 
  # scale_alpha_manual(values = c(1,1), na.value = .5) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_blank(), 
        axis.ticks = element_blank()) + 
  xlab("") + ylab("")
