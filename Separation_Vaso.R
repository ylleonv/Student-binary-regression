library(tidyr);library(dplyr)
library(ggplot2)

set.seed(2)

x1 <- rnorm(11,1,1)
x2 <- rnorm(11,0,1)

dat1 <- tbl_df(cbind(x1,x2))
# dat1[12,] <- t(c(1,0))
# dat1[13,] <- t(c(2,0.8))

ggplot(dat1, aes(x=x1,y=x2))+
  geom_text(
    label=rownames(dat1),
    check_overlap = T
  )+
  geom_abline(aes(intercept = -0.8, slope = 0.8))

group1 <- c(1,4,6,8,10)
dat1[group1,"y"] <- 1
dat1[rownames(dat1)[-group1],"y"] <- 0

ggplot(dat1, aes(x=x1,y=x2, color = factor(y)))+
  geom_text(
    label=rownames(dat1),
    check_overlap = T
  )+
  geom_abline(aes(intercept = -0.8, slope = 0.8))

glm(y~x1+x2,dat1,family = binomial(link = "logit"))

dat1[11,"y"] <- 1

ggplot(dat1, aes(x=x1,y=x2, color = factor(y)))+
  geom_text(
    label=rownames(dat1),
    check_overlap = T
  )+
  geom_abline(aes(intercept = -0.8, slope = 0.8))

mod2 <- glm(y~x1+x2,dat1,family = binomial(link = "logit"))

ux <- predict(mod2,type = "link")

dat2 <- tbl_df(cbind(ux,y=dat1$y))
ggplot(data = dat2, aes(x=ux,y=y))+
  geom_text(
    label=rownames(dat1),
    check_overlap = T
  )















# DATA VASO ---------------------------------------------------------------

data(vaso, package="catdata")
vaso <- tbl_df(vaso)
vaso$vaso[vaso$vaso==2] <- 0

plot1 <- ggplot(data = vaso, aes(x = vol, y = rate, color=factor(vaso))) +
  geom_point()+
  geom_text(
    data = vaso[c(4,18,29),], aes(x = vol, y = rate),
    label = c(4,18,29), hjust = 0.5, vjust = 1.3
    # check_overlap = T
  )+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none")+
  xlab("ln of Volume") + ylab("ln of Rate")+
  geom_abline(intercept = 0.77, slope = -1.1, linetype="dashed")+
  scale_color_manual("", labels = c("0", "1"), values = c("#E69F00", "#00AFBB"))

plot1
# tikz(file = "separation_vaso.tex", standAlone=F,width = 6.5, height = 3)
# endoffile <- dev.off() 

vaso[c(4,18,29),"vaso"] = 1
# vaso[c(29),"vaso"] = 1
mod0 <- glm(factor(vaso)~vol+rate, data = vaso,
            family = binomial(link = "logit"))
logLik(mod0)
ux <- predict(mod0,type = "link")
# ux[4] <- -3.1

ux1 <- predict(mod0,type = "response")

vaso1 <- vaso
vaso1[c(4,18,29),"vaso"] = 0
mod1 <- glm(factor(vaso)~vol+rate, data = vaso1,
            family = binomial(link = "logit"))
# ux1 <- predict(mod1,type = "link")
ux11 <- predict(mod1,type = "response")

dat2 <- tbl_df(cbind(ux,y=vaso$vaso))
dat2$color <- 1
dat2$color[c(4,18,29)] <- 3

dat3 <- dat2[c(4,18,29),]
dat3$tag <- c(4,18,29)
dat3$y <- c(0,0,0)
dat3[1,1] <- -3.1

dat4 <- dat2[c(4,18,29),]
dat4[1,1] <- -3.1

plot2 <- ggplot(data = dat2, aes(x=ux,y=y))+
  geom_point(data=dat4,aes(x=ux,y=y),color="#00AFBB")+
  geom_point(data = dat2[-c(4,18,29),],aes(x=ux,y=y),color="gray")+
  geom_point(data= dat3, aes(x=ux,y=y), color="#E69F00")+
  geom_text(
    # size = 3,
    color="#E69F00",
    data = dat3, aes(x = ux, y = y),
    label = c(4,18,29), hjust = 0.5, vjust = 1.3
    # check_overlap = T
  )+
  geom_text(
    # size = 2,
    color="#00AFBB",
    data = dat4, aes(x = ux, y = y,color=factor(color)),
    label = c(4,18,29), hjust = 0.5, vjust = 1.3
    # check_overlap = T
  )+
  geom_line(aes(x=ux,y=ux1),color="#00AFBB")+
  geom_line(aes(x=ux,y=ux11),color="#E69F00")+
  # geom_vline(xintercept = 0.65)+
  # scale_color_manual("", labels = c("0", "1"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none")+
  xlab("xu") + ylab("y") + ylim(-0.03,1)
plot2
# tikz(file = "vaso_sep.tex", standAlone=F,width = 6.5, height = 3)
# endoffile <- dev.off() 

library(ggpubr)
ggarrange(plot1,plot2,ncol = 2, nrow = 1, labels = c("",""))

