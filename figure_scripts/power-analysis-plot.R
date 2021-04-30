setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(foreach)
library(gridExtra)
library(reshape2)
library(EnvStats)
# source("resampling-test.R")
theme_set(theme_bw())

ts.frame <- read.table("resampling_data.txt", header = TRUE) # data frame containing test statistics (both from actual data and from resampled data)
ts.lscm.data <- subset(ts.frame, b==0 & model == "lscm") # LSCM test statistic on actual data
ts.noconf.data <- subset(ts.frame, b==0 & model == "noconf") # difference in sample avg on actual data
# ts.lscm.data <- -0.018
# ts.noconf.data <- 0.073

power.frame <- read.table("2021-02-15_power-analisys.txt")

power.frame.agg <- aggregate(pv ~ b + test, FUN = function(pv) mean(pv<.05), power.frame)
power.frame.agg$test <- factor(power.frame.agg$test, 
                       levels=c("t-test", "resampling-conconf", "resampling-adjustH", "resampling-lscm"), 
                       labels = c("two-sample t-test", "resampling - no adjustment", "resampling - adjust for W", "resampling - LSCM"))

q95 <- qbinom(p=.95, size=100, prob=.05)
p.power <- ggplot(power.frame.agg, aes(b,100*pv,col=test, lty = test)) + 
  geom_rect(xmin=-0.02,xmax=.1,ymin = 0, ymax=q95, col="transparent", fill="grey", alpha=.03) + 
  geom_line(size=1) + 
  scale_color_manual(values = c("#CC79A7", "#E69F00", "#0072B2", "#009E73")) + 
  scale_linetype_manual(values = c(1,2,1,1)) + 
  # geom_hline(yintercept=5, lty=2, col="red") + 
  geom_vline(xintercept = abs(ts.lscm.data), lty = "dashed") + 
  geom_vline(xintercept = abs(ts.noconf.data), lty = "dotted") +
  xlab("causal effect of X (beta)") + 
  # ylab("empirical statistical power") + 
  ylab(expression(paste("number of rejections of ", H[0]))) + 
  theme(legend.position = c(0.65,0.35))
p.power

pdf("power.pdf", width = 6, height = 4)
p.power
dev.off()

pdf("../figures/power.pdf", width = 6, height = 4)
p.power
dev.off()

