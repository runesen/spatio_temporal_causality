setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(gridExtra)
theme_set(theme_bw())

## H

set.seed(1)
n <- 47
hsim <- as.numeric(arima.sim(model=list(ar=.9,sd=.5),n=n))
plot(hsim,type="l")
pt <- numeric(n)
pt[c(10,24,38)] <- 1
df <- data.frame(x = 1:n, 
                 h = hsim, 
                 points = pt)
plot(1:n,hsim,type="l")

pH <- ggplot(df, aes(x,h)) + geom_line(lty="dashed") + geom_point(aes(alpha=I(pt)),size=2) + xlab("") + ylab("") + 
  theme(panel.grid = element_blank(), 
        axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(), 
        plot.margin = unit(c(0,-6,0,-12),"mm"))

## X and Y

set.seed(2)

nn <- 11
m <- 100
h <- rep(seq(0,2,length.out=nn),each=m)
x <- h + .5*rnorm(nn*m)
y <- 2*h+.25*h*x + .3*rnorm(nn*m)
df1 <- data.frame(x=x,y=y,h=h)

bconf <- coefficients(lm(y~x))
b05 <- coefficients(lm(y~x, data=subset(df1, abs(h-.4)<.05)))
b10 <- coefficients(lm(y~x, data=subset(df1,h==1)))
b15 <- coefficients(lm(y~x, data=subset(df1,abs(h-1.6)<.05)))
bavgcausal <- colMeans(matrix(coefficients(lm(y~x*factor(h)))[c(1,3:12,2,13:22)],nrow=nn))
btrue <- c(2, .25) # ground truth if h is uniform(0,2)

pconf <- ggplot(df1, aes(x,y)) + geom_point(alpha=.7) + 
  xlab("X") + ylab("Y") + 
  coord_cartesian(xlim=range(x),ylim=range(y)) + 
  geom_abline(intercept = bconf[1], slope = bconf[2], col="red", size=1.5) + 
  geom_abline(intercept = btrue[1], slope = btrue[2], col="green", size=1.5, lty = 2) + 
  theme(panel.grid = element_blank(), 
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_text(size=15))

pXY1 <- ggplot(subset(df1, abs(h-1) < .05), aes(x,y)) + geom_point(alpha=.7) + 
  # xlab("X") + ylab("Y") + 
  coord_cartesian(xlim=range(x),ylim=range(y)) + 
  geom_abline(intercept = b10[1], slope = b10[2], col="blue", size=1) + 
  theme(panel.grid = element_blank(), 
        axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())

pXY2 <- ggplot(subset(df1, abs(h-.4)<.05), aes(x,y)) + geom_point(alpha=.7) + 
  # xlab("X") + ylab("Y") + 
  coord_cartesian(xlim=range(x),ylim=range(y)) + 
  geom_abline(intercept = b05[1], slope = b05[2], col="blue", size=1) + 
  theme(panel.grid = element_blank(), 
        axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())

pXY3 <- ggplot(subset(df1, abs(h-1.6)<.05), aes(x,y)) + geom_point(alpha=.7) + 
  # xlab("X") + ylab("Y") + 
  coord_cartesian(xlim=range(x),ylim=range(y)) + 
  geom_abline(intercept = b15[1], slope = b15[2], col="blue", size=1) + 
  theme(panel.grid = element_blank(), 
        axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())

pavgcausal <- ggplot(df1, aes(x,y)) + geom_point(alpha=.7) + 
  xlab("X") + ylab("Y") + 
  coord_cartesian(xlim=range(x),ylim=range(y)) + 
  geom_abline(intercept=bavgcausal[1], slope=bavgcausal[2], col="blue", size=1.5) + 
  geom_abline(intercept = btrue[1], slope = btrue[2], col="green", size=1.5, lty = 2) + 
  theme(panel.grid = element_blank(), 
        axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_text(size=15))

grid.arrange(pconf, pXY1, pXY2, pXY3, pavgcausal, nrow=1)

pdf("../figures/pH_new.pdf", width = 10, height = 1.5)
pH
dev.off()

pdf("../figures/pconf.pdf", width = 4, height = 3)
pconf
dev.off()

pdf("../figures/pavgcausal.pdf", width = 4, height = 3)
pavgcausal
dev.off()


pdf("../figures/pXY1_new.pdf", width = 4, height = 3)
pXY1
dev.off()

pdf("../figures/pXY2_new.pdf", width = 4, height = 3)
pXY2
dev.off()

pdf("../figures/pXY3_new.pdf", width = 4, height = 3)
pXY3
dev.off()
