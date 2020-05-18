setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(foreach)
library(ggplot2)
theme_set(theme_bw())

m <- 511
x <- seq(1,m,.5)
y <- floor(log2(x)) %% 2 + 1
df <- data.frame(x=x,y=y,z=floor(log2(x)))
p <- ggplot(df, aes(x,y,group=z)) + 
  geom_line() + 
  #geom_point(size=.1) + 
  xlab("t") + 
  # ylab(expression(paste("pairity(floor(",log[2](t),"))"))) + 
  ylab(expression(paste("parity [",log[2](t),"]"))) + 
  scale_x_continuous(breaks = 2^(2:9)) + 
  scale_y_continuous(breaks = c(1,2), labels = c("even", "odd")) + 
  theme(axis.text = element_text(size=10))
p

pdf("../figures/plot_remark.pdf", width=15, height=2)
p
dev.off()