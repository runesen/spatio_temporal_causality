setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# dat_meta <- read.table("uniqueID_country-province.txt", header = TRUE, sep="\t")
dat_meta <- read.csv("uniqueID_country_Biome_LS_Anthrome_20190130.csv", header=T, sep= ",", dec=".")
dat_popdens <- read.csv("GPWv4_summary_20190129.csv") 
dat_FL <- read.csv("Hansen_Summaries_ALL_20190821.csv")
dat_conflict <- read.csv("_PRIO-summaries_ALL_bestEstimate_20191219.csv", header = T, sep=",", dec=".")
dat_roaddist <- read.csv("RoadDist_Summaries.csv", header = TRUE)

head(dat_popdens)
head(dat_FL)
head(dat_conflict, 22)
head(dat_meta)
head(dat_roaddist)

# aggregating onto common grid
shared_polygons <- Reduce(intersect, list(unique(dat_meta$UniqueID), 
                                          unique(dat_popdens$UniqueID), 
                                          unique(dat_FL$UniqueID), 
                                          unique(dat_conflict$PolygonID), 
                                          unique(dat_roaddist$UniqueID))
                          )

dat_meta <- subset(dat_meta, UniqueID %in% shared_polygons)
dat_popdens <- subset(dat_popdens, UniqueID %in% shared_polygons)
dat_FL <- subset(dat_FL, UniqueID %in% shared_polygons)
dat_conflict <- subset(dat_conflict, PolygonID %in% shared_polygons)
dat_roaddist <- subset(dat_roaddist, UniqueID %in% shared_polygons)

# output data frame
dat <- dat_conflict[,c("PolygonID", "Year", "nr_fat_all")]
colnames(dat)[3] <- "nr_fatalities"

nt <- length(unique(dat$Year))
ns <- length(unique(dat$PolygonID))
n <- ns*nt

# sorting
dat <- dat[order(dat$PolygonID, dat$Year),]
dat_meta <- dat_meta[order(dat_meta$UniqueID),]
dat_FL <- dat_FL[order(dat_FL$UniqueID),]
dat_conflict <- dat_conflict[order(dat_conflict$PolygonID, dat_conflict$Year),]
dat_roaddist <- dat_roaddist[order(dat_roaddist$UniqueID),]

# year 2000 = forest loss NA
FL_matrix <- cbind(rep(NA, ns), as.matrix(dat_FL[,3:(2+nt-1)]))
dat$FL_km <- c(t(FL_matrix))
# population density data updated 2000, 2005, 2010, 2015
popdens_matrix <- as.matrix(dat_popdens[,rep(2:5, c(5,5,5,4))])
dat$PopDens <- c(t(popdens_matrix))
# forest in year 2000
dat$xFor2000_ge25 <- rep(dat_FL$F2000_km_th25, each = nt)
# distance to road
dat$RoadDist <- rep(dat_roaddist$RoadDist, each = nt)

# latitude, longituge
dat$lon <- rep(dat_meta$POINT_X, each = nt)
dat$lat <- rep(dat_meta$POINT_Y, each = nt)
dat$country_name <- rep(dat_meta$GID_0, each = nt)


## Cumulative forest loss
l <- lapply(1:(nt-1), function(i){
  v <- c(rep(0,i),dat$FL_km[-c((n-i+1):n)])
  ind <- c(sapply(1:i, function(j) seq(j, n-(nt-j), by=nt)))
  v[ind] <- 0
  v
})
A <- cbind(dat$FL_km, Reduce(cbind,l))
FL_km_cum <- rowSums(A,na.rm=1)


## Forest cover
dat$xFor_ge25 <- dat$xFor2000_ge25 - FL_km_cum
dat$xFor_ge25 <- pmax(0, dat$xFor_ge25)

## Set all FL_km to NA for which there was no forest left in the year before
forest_zero <- which(dat$xFor_ge25==0)
set_FL_to_NA <- 1+setdiff(forest_zero, seq(nt,n,by=nt))
dat$FL_km[set_FL_to_NA] <- NA


dat[dat$country_name == "",]$country_name <- NA

dat0 <- subset(dat, country_name == "COL")
write.table(dat0, "data_xy_colombia_20200616.txt", quote = FALSE)
