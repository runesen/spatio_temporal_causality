# construct summary data set
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(gridExtra)
library(viridis)
theme_set(theme_bw())

df <- read.table("data_xy_colombia_20200616.txt", header = TRUE)
df.noforest <- subset(df, xFor2000_ge25 == 0 & Year == 2000)
poly.noforest <- df.noforest$PolygonID
df.agg <- aggregate(FL_km ~ lon + lat, df, sum, na.rm = TRUE, na.action=NULL)
colnames(df.agg)[3] <- "FL_km_cumulative"
df.agg$PolygonID <- aggregate(PolygonID ~ lon + lat, df, function(x) x[1])$PolygonID
df.agg$FL_km_cumulative[df.agg$PolygonID %in% poly.noforest] <- NA
df.agg$popdens_average <- aggregate(PopDens ~ lon + lat, df, mean)$PopDens
df.agg$RoadDist <- aggregate(RoadDist ~ lon + lat, df, mean)$RoadDist
df.agg$conflict_location <- aggregate(nr_fatalities ~ lon + lat, df, function(x) any(x>0))$nr_fatalities

c <- df$nr_fatalities>0
fl <- df$FL_km
mean(fl[c], na.rm=1) # 0.266
mean(fl[!c], na.rm=1) # 0.193
round(t.test(fl~c)$p.value,12) # 8.89e-10

write.table(df.agg, "data_xy_colombia_temporally_aggregated_20200616.txt", quote=FALSE)


################# summaries on province level

df.prov <- read.table("uniqueID_country-province.txt", header = TRUE, sep = "\t")
df.prov <- subset(df.prov, GID_0 == "COL")
ns <- length(unique(df$PolygonID))
nt <- length(unique(df$Year))
df <- df[order(df$PolygonID),]
df.prov <- df.prov[order(df.prov$UniqueID),]
df$province <- factor(rep(as.character(df.prov$NAME_1), each = nt))

df.agg.prov <- aggregate(FL_km ~ province, df, sum, na.rm = TRUE, na.action=NULL)
df.agg.prov$nr_conflicts <- aggregate(nr_fatalities ~ province, df, function(x) sum(x>0))$nr_fatalities

write.table(df.agg.prov, "data_xy_colombia_regional_summaries_20200432.txt", quote=FALSE)

