
# Influence function plot -------------------------------------------------

eta <- seq(-5,5,0.01)
logis <- plogis(eta)
normal <- pnorm(eta)
cauchy <- pcauchy(eta)
student <- pt(eta,df = 0.3)
# gompertz <-  1 - exp(-(exp(eta)- 1))
gompertz <-  1 - exp(-(exp(eta)))
gumbel <- exp(-exp(-eta))
# gompertz <- flexsurv::pgompertz(eta, shape = 1, rate = 1, lower.tail = TRUE, log.p = FALSE)

links <- rbind(data.frame(eta,value=logis,Link="Logistic"),
               data.frame(eta,value=normal,Link="Normal"),
               data.frame(eta,value=cauchy,Link="Cauchy"),
               # data.frame(eta,value=student,Link="Student(0.3)"),
               data.frame(eta,value=gompertz,Link="Gompertz"),
               data.frame(eta,value=gumbel,Link="Gumbel"))

summary(links)

ggplot(data = links)+
  geom_line(aes(x = eta, y = value, color=factor(Link))) +
  scale_color_manual(values = c("#e6194B", "#3cb44b","#f58231","#4363d8","#42d4f4"))+
  labs(title="", x="$eta$", y="$pi$", color = "Link function")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),) 

tikz(file = "tex/State of the art/Link functions.tex", standAlone=F,width = 5, height = 3)

endoffile <- dev.off() 


# Normalization -----------------------------------------------------------

normalization <- 0.95
qp <- qlogis(normalization)
qgompertz <- function(p){log(1-log(1-p))}
qgumbel <- function(p){-log(-log(p))}

s0_norm <- qp/(qnorm(normalization)-qnorm(0.5))
s0_cauchy <- qp/(qcauchy(normalization)-qcauchy(0.5))
s0_t <- qp/(qt(normalization,0.3)-qt(0.5,0.3))
s0_gompertz <- qp/(qgompertz(normalization)-qgompertz(0.5))
s0_gumbel <- qp/(qgumbel(normalization)-qgumbel(0.5))
m0_gompertz <- (qgompertz(0.5)*qp)/(qgompertz(normalization)-qgompertz(0.5))
m0_gumbel <- (qgumbel(0.5)*qp)/(qgumbel(normalization)-qgumbel(0.5))

logis <- plogis(eta)
normal <- pnorm(eta/s0_norm)
cauchy <- pcauchy(eta/s0_cauchy)
student <- pt(eta/s0_t,df = 0.3)
gompertz <- (1-exp(-exp((eta+m0_gompertz)/s0_gompertz)))
gumbel <- (exp(-exp(-((eta+m0_gumbel)/s0_gumbel))))

links <- rbind(data.frame(eta,value=logis,Link="Logistic"),
               data.frame(eta,value=normal,Link="Normal"),
               data.frame(eta,value=cauchy,Link="Cauchy"),
               data.frame(eta,value=student,Link="Student(0.3)"),
               # data.frame(eta,value=gompertz,Link="Gompertz"),
               data.frame(eta,value=gumbel,Link="Gumbel"))


plot095 <- ggplot(data = links)+
  geom_line(aes(x = eta, y = value, color=factor(Link))) +
  labs(title="", x="eta", y="pi", color = "Link function")+
  scale_color_manual(values = c("#e6194B", "#3cb44b","#f58231","#4363d8","#42d4f4"))+
  geom_segment(aes(x = 2.94, xend=2.94, y=-0.1,yend=0.95), linetype = "dashed", color = "#a9a9a9")+
  geom_segment(aes(x = 0, xend=0, y=-0.1,yend=0.5), linetype = "dashed", color = "#a9a9a9")+
  geom_segment(aes(x = -6.9, xend=2.94, y=0.95,yend=0.95), linetype = "dashed", color = "#a9a9a9")+
  geom_segment(aes(x = -6.9, xend=0, y=0.5,yend=0.5), linetype = "dashed", color = "#a9a9a9")+
  scale_y_continuous("", c(0,0.5,0.95), labels=c("0","0.5","0.95"))+
  scale_x_continuous("", c(0,2.94), labels=c("0","2.94"))+
  coord_cartesian(ylim=c(0,1),xlim=c(-4.5,4.5))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none") 


normalization <- 0.73
qp <- qlogis(normalization)
s0_norm <- qp/(qnorm(normalization)-qnorm(0.5))
s0_cauchy <- qp/(qcauchy(normalization)-qcauchy(0.5))
s0_t <- qp/(qt(normalization,0.3)-qt(0.5,0.3))
s0_gompertz <- qp/(qgompertz(normalization)-qgompertz(0.5))
s0_gumbel <- qp/(qgumbel(normalization)-qgumbel(0.5))
m0_gompertz <- (qgompertz(0.5)*qp)/(qgompertz(normalization)-qgompertz(0.5))
m0_gumbel <- (qgumbel(0.5)*qp)/(qgumbel(normalization)-qgumbel(0.5))

logis <- plogis(eta)
normal <- pnorm(eta/s0_norm)
cauchy <- pcauchy(eta/s0_cauchy)
student <- pt(eta/s0_t,df = 0.3)
gompertz <- (1-exp(-exp((eta+m0_gompertz)/s0_gompertz)))
gumbel <- (exp(-exp(-((eta+m0_gumbel)/s0_gumbel))))

links <- rbind(data.frame(eta,value=logis,Link="Logistic"),
               data.frame(eta,value=normal,Link="Normal"),
               data.frame(eta,value=cauchy,Link="Cauchy"),
               data.frame(eta,value=student,Link="Student(0.3)"),
               # data.frame(eta,value=gompertz,Link="Gompertz"),
               data.frame(eta,value=gumbel,Link="Gumbel"))


plot073 <- ggplot(data = links)+
  geom_line(aes(x = eta, y = value, color=factor(Link))) +
  labs(title="", x="eta", y="pi", color = "Link function")+
  scale_color_manual(values = c("#e6194B", "#3cb44b","#f58231","#4363d8","#42d4f4"))+
  geom_segment(aes(x = 1, xend=1, y=-0.1,yend=0.73), linetype = "dashed", color = "#a9a9a9")+
  geom_segment(aes(x = 0, xend=0, y=-0.1,yend=0.5), linetype = "dashed", color = "#a9a9a9")+
  geom_segment(aes(x = -6.9, xend=qp, y=0.73,yend=0.73), linetype = "dashed", color = "#a9a9a9")+
  geom_segment(aes(x = -6.9, xend=0, y=0.5,yend=0.5), linetype = "dashed", color = "#a9a9a9")+
  scale_y_continuous("", c(0,0.5,0.73), labels=c("0","0.5","0.73"))+
  scale_x_continuous("", c(0,1), labels=c("0","1"))+
  coord_cartesian(ylim=c(0,1),xlim=c(-4.5,4.5))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none") 


normalization <- 0.99
qp <- qlogis(normalization)
s0_norm <- qp/(qnorm(normalization)-qnorm(0.5))
s0_cauchy <- qp/(qcauchy(normalization)-qcauchy(0.5))
s0_t <- qp/(qt(normalization,0.3)-qt(0.5,0.3))
s0_gompertz <- qp/(qgompertz(normalization)-qgompertz(0.5))
s0_gumbel <- qp/(qgumbel(normalization)-qgumbel(0.5))
m0_gompertz <- (qgompertz(0.5)*qp)/(qgompertz(normalization)-qgompertz(0.5))
m0_gumbel <- (qgumbel(0.5)*qp)/(qgumbel(normalization)-qgumbel(0.5))

logis <- plogis(eta)
normal <- pnorm(eta/s0_norm)
cauchy <- pcauchy(eta/s0_cauchy)
student <- pt(eta/s0_t,df = 0.3)
gompertz <- (1-exp(-exp((eta+m0_gompertz)/s0_gompertz)))
gumbel <- (exp(-exp(-((eta+m0_gumbel)/s0_gumbel))))

links <- rbind(data.frame(eta,value=logis,Link="Logistic"),
               data.frame(eta,value=normal,Link="Normal"),
               data.frame(eta,value=cauchy,Link="Cauchy"),
               data.frame(eta,value=student,Link="Student(0.3)"),
               # data.frame(eta,value=gompertz,Link="Gompertz"),
               data.frame(eta,value=gumbel,Link="Gumbel"))

plot099 <- ggplot(data = links)+
  geom_line(aes(x = eta, y = value, color=factor(Link))) +
  labs(title="", x="eta", y="pi", color = "Link function")+
  scale_color_manual(values = c("#e6194B", "#3cb44b","#f58231","#4363d8","#42d4f4"))+
  geom_segment(aes(x = 4.6, xend=4.6, y=-0.1,yend=0.99), linetype = "dashed", color = "#a9a9a9")+
  geom_segment(aes(x = 0, xend=0, y=-0.1,yend=0.5), linetype = "dashed", color = "#a9a9a9")+
  geom_segment(aes(x = -6.9, xend=4.6, y=0.99,yend=0.99), linetype = "dashed", color = "#a9a9a9")+
  geom_segment(aes(x = -6.9, xend=0, y=0.5,yend=0.5), linetype = "dashed", color = "#a9a9a9")+
  scale_y_continuous("", c(0,0.5,0.99), labels=c("0","0.5","0.99"))+
  scale_x_continuous("", c(0,4.6), labels=c("0","4.6"))+
  coord_cartesian(ylim=c(0,1),xlim=c(-4.5,4.5))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none") 


library(ggpubr)

ggarrange(plot073,plot095,plot099, ncol = 3, nrow = 1, labels = c("","",""))

tikz(file = "normalization_shapes.tex", standAlone=F,width = 7, height = 3)
endoffile <- dev.off() 
