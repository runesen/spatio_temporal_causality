###################################
# computing p-values for permutation 
# tests and plots 
##################################"
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(foreach)
library(gridExtra)
library(reshape2)
theme_set(theme_bw())


res.frame <- read.table("resampling_data.txt", header = TRUE)
subset(res.frame, b == 0)
# ts_noconf = 0.073
# ts_popdens = 0.037
# ts_roaddist = 0.037
# ts_lscm = -0.018

# p-value for two-sided test
pvfun <- function(ts, na.rm){
  ts.data <- ts[1]
  ts.res <- ts[-1]
  ts.res <- ts.res[!is.na(ts.res)]
  B <- length(ts.res)
  min(1,2*min((1+sum(ts.res <= ts.data))/(1+B), (1+sum(ts.res >= ts.data))/(1+B)))
} 
pvframe <- aggregate(ts ~ model, res.frame, FUN = pvfun)
colnames(pvframe)[2] <- "pval"


pdf("../figures/teststat_resampled_noconf.pdf", width=5, height = 3.5)
ggplot(subset(res.frame, b > 0 & model == "noconf"), aes(x=ts)) + geom_histogram(color="gray", alpha=0.7) + 
  geom_vline(xintercept=subset(res.frame, model == "noconf")$ts[1], col="red", size = 1) + 
  ggtitle(paste0("p-value = ", sprintf("%.3f", subset(pvframe, model == "noconf")$pval))) + 
  xlab("test statistic") + theme(text = element_text(size=15), plot.title = element_text(hjust=.5))
dev.off()

pdf("../figures/teststat_resampled_popdens.pdf", width=5, height = 3.5)
ggplot(subset(res.frame, b > 0 & model == "popdens"), aes(x=ts)) + geom_histogram(color="gray", alpha=0.7) + 
  geom_vline(xintercept=subset(res.frame,  model == "popdens")$ts[1], col="red", size = 1) + 
  ggtitle(paste0("p-value = ", sprintf("%.3f", subset(pvframe,  model == "popdens")$pval))) + 
  xlab("test statistic") + theme(text = element_text(size=15), plot.title = element_text(hjust=.5))
dev.off()

pdf("../figures/teststat_resampled_roaddist.pdf", width=5, height = 3.5)
ggplot(subset(res.frame, b > 0 & model == "roaddist"), aes(x=ts)) + geom_histogram(color="gray", alpha=0.7) + 
  geom_vline(xintercept=subset(res.frame, model == "roaddist")$ts[1], col="red", size = 1) + 
  ggtitle(paste0("p-value = ", sprintf("%.3f", subset(pvframe, model == "roaddist")$pval))) + 
  xlab("test statistic") + theme(text = element_text(size=15), plot.title = element_text(hjust=.5))
dev.off()

pdf("../figures/teststat_resampled_lscm.pdf", width=5, height = 3.5)
ggplot(subset(res.frame, b > 0 & model == "lscm"), aes(x=ts)) + geom_histogram(color="gray", alpha=0.7) + 
  geom_vline(xintercept=subset(res.frame, model == "lscm")$ts[1], col="red", size = 1) + 
  ggtitle(paste0("p-value = ", sprintf("%.3f", subset(pvframe, model == "lscm")$pval))) + 
  xlab("test statistic") + theme(text = element_text(size=15), plot.title = element_text(hjust=.5))
dev.off()

