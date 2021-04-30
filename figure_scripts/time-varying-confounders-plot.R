setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(gridExtra)
library(fields)
library(foreach)
library(dplyr, warn.conflicts = FALSE)
theme_set(theme_bw())


df.coef <- read.table("time-var-conf.txt", header = TRUE)
df.coef$b0[df.coef$var==0 & df.coef$method == "adj-for-all-H"] <- NA
df.coef$method <- factor(df.coef$method, 
                         levels = c("adj-for-time-inv-H", "adj-for-all-H", "marginal"), 
                         labels = c("LSCM", "LSCM~+~Hbar", "no~adjustm."))
df.coef$var <- factor(df.coef$var, 
                      levels = c(0,0.05,0.1,0.15,0.2),
                      labels = c(expression(paste(sigma^2 , " = 0.00")),
                                 expression(paste(sigma^2 , " = 0.05")),
                                 expression(paste(sigma^2 , " = 0.10")),
                                 expression(paste(sigma^2 , " = 0.15")),
                                 expression(paste(sigma^2 , " = 0.20"))))


pdf("../figures/b0hat.pdf", width = 10, height = 3.5)
ggplot(df.coef, aes(b0)) + geom_histogram() + 
  facet_grid(method~var, scales = "free_y", labeller = label_parsed) + 
  geom_vline(xintercept=1, col="red") + 
  xlab(expression(paste("estimate of ", beta[0]))) 
dev.off()

pdf("../figures/b1hat.pdf", width = 10, height = 3.5)
ggplot(df.coef, aes(b1)) + geom_histogram() + 
  facet_grid(method~var, scales = "free_y", labeller = label_parsed) + 
  geom_vline(xintercept=2, col="red") + 
  xlab(expression(paste("estimate of ", beta[1])))
dev.off()


