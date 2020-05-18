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

# dat <- read.table("../Conflict_Data/data_xy_20190205.txt", header = TRUE)
# dat1 <- subset(dat, country_name == "COL")
# df <- read.table("data_xy_colombia_20191219.txt", header = T)
df <- read.table("data_xy_colombia_20200327.txt", header = T)
df.a <- aggregate(nr_fatalities ~ Year, df, sum)
df.a$FL_km <- c(NA, aggregate(FL_km ~ Year, subset(df, Year %in% 2001:2018), sum, na.rm=1)$FL_km)

rects <- data.frame(xstart = c(2000, 2012, 2016), xend = c(2012, 2016, 2018), 
                    col = factor(c(1,2,1)))

p.c <- ggplot() + geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=col), alpha=.3) + 
  scale_fill_manual(values = c("transparent", "red")) + 
  geom_point(data = df.a, aes(Year, nr_fatalities)) + 
  geom_line(data = df.a, aes(Year, nr_fatalities), size = .8) + 
  geom_vline(xintercept = 2016, col = "red", size=1.2) + 
  xlab("year") + ylab("number of conflicts") + 
  theme(legend.position = "none")
p.c

p.fl <- ggplot() + geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=col), alpha=.3) + 
  scale_fill_manual(values = c("transparent", "red")) + 
  geom_point(data=df.a, aes(Year-.5, FL_km)) + 
  geom_line(data=df.a, aes(Year-.5, FL_km), size=.8) + 
  geom_vline(xintercept = 2016, col = "red", size=1.2) + 
  xlab("year") + ylab("forest loss (km2)") + 
  theme(legend.position = "none")
p.fl

pdf("../figures/effect_of_peace_agreement1.pdf", width=6, height = 4)
p.c
dev.off()

pdf("../figures/effect_of_peace_agreement2.pdf", width=6, height = 4)
p.fl
dev.off()

