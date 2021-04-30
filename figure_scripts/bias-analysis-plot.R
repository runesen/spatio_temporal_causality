setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(foreach)
library(gridExtra)
library(reshape2)
library(EnvStats)
theme_set(theme_bw())

out <- read.table("2021-02-19_bias-analisys.txt", header = TRUE)
out$method <- factor(out$method, levels = c("conf", "adjustH", "lscm"))
out.agg <- aggregate(bias ~ a+b+method, FUN=mean, out)
out.agg$ts.true <- aggregate(ts.true ~ a+b+method, FUN=mean, out)$ts.true
out.agg$method <- factor(out.agg$method, 
                         levels = c("conf", "adjustH", "lscm"),
                         labels = c("no adjustment", "adjustment for W", "LSCM"))


p.bias <- ggplot(out.agg, aes(a,b,fill=bias)) + 
  geom_tile() + facet_grid(.~method) + 
  scale_fill_gradient2() + 
  geom_text(aes(a,b,label=round(bias,3))) + 
  xlab("causal effect of log W (alpha)") + ylab("causal effect of X (beta)")
p.bias

pdf("bias.pdf", width = 15*.85, height = 4.9*.85)
p.bias
dev.off()

pdf("../figures/bias.pdf", width = 15*.85, height = 4.9*.85)
p.bias
dev.off()