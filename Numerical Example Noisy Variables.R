# Numerical separation ----------------------------------------------------

S<-function(x,s){ return  (1/(1+exp(-x/s)))}

# Illustrative numerical example

set.seed(23)
n_obs <- 100

x <- rnorm(n_obs)
x

alpha <- 0
beta1 <- 1

eta1 <- alpha + as.matrix(x) %*% (beta1)
s <-  0.5
h_x_beta <- S(eta1,s)
y <- rbinom(n = n_obs,size = 1,prob = h_x_beta)  
y[which.min(h_x_beta)]=1
dat <- data.frame(y,x)
dat$y <- as.factor(dat$y)
plot(dat$x,dat$y)

mod_log <- glm(formula = y ~ x , data = dat, family = binomial(link = "logit"))
eta_pred <- predict(mod_log)
dat_pred <- data.frame(y,eta_pred)
plot(dat_pred$eta_pred,dat_pred$y)
plot(dat$x,predict(mod_log, type = "response"))

coef_log <- coef(mod_log)
ggplot(dat,aes(x=x,y=0, shape=y,color=y))+
  geom_point()+
  geom_vline(xintercept = -coef_log[1]/coef_log[2])


mod_stu <- glm(formula = y ~ x, data = dat, family = binomial(link = Gosset(0.8)))
eta_pred_stu <- predict(mod_stu)
dat_pred_stu <- data.frame(y,eta_pred_stu)
plot(dat_pred_stu$eta_pred_stu,dat_pred_stu$y)
plot(dat$x,predict(mod_stu, type = "response"))

coef_stu <- coef(mod_stu)
ggplot(dat,aes(x=x,y=0, shape=y,color=y))+
  geom_point()+
  geom_vline(xintercept = -coef_stu[1]/coef_stu[2])

set.seed(25)
dat2 <- data.frame(y,x,n=rnorm(n_obs))
dat2$y <- as.factor(dat2$y)

mod_log_n1 <- glm(formula = y ~ x + n, data = dat2, family = binomial(link = "logit"))
eta_pred <- predict(mod_log_n1)
dat2_pred <- data.frame(y,eta_pred)
plot(dat2_pred$eta_pred,dat2_pred$y)
plot(dat2$x,predict(mod_log_n1, type = "response"))

coef_logn1 <- coef(mod_log_n1)
ggplot(dat2,aes(x=x,y=n, shape=y,color=y))+
  geom_point()+
  geom_segment(aes(x = 0, y = -coef_logn1[1]/coef_logn1[3],
    xend = -coef_logn1[1]/coef_logn1[2], yend = 0))+
  ylim(-4,3)
# +
#   geom_abline(slope = 27.6428, intercept = -2.78055)

mod_stu_n1 <- glm(formula = y ~ x + n , data = dat2, family = binomial(link = Gosset(0.8)))
eta_pred_stu_n1 <- predict(mod_stu_n1)
dat2_pred_stu <- data.frame(y,eta_pred_stu_n1)
plot(dat2_pred_stu$eta_pred_stu_n1,dat2_pred_stu$y)
plot(dat2$x,predict(mod_stu_n1, type = "response"))

coef_stun1 <- coef(mod_stu_n1)
ggplot(dat2,aes(x=x,y=n, shape=y,color=y))+
  geom_point()+
  geom_segment(aes(x = 0, y = -coef_stun1[1]/coef_stun1[3],
    xend = -coef_stun1[1]/coef_stun1[2], yend = 0))+
  ylim(-4,3)
# +
#   geom_abline(slope = 25.2964, intercept = -2.449327)



# Numerical illustration --------------------------------------------------

plot1_po_pn <- plot1_po_n <- plot1_o_pn <- plot1_o_n <- NULL
seedd <- seedd+1
seedd <- 50
# 50 71
set.seed(seedd*8) # 4  ou 11
{
  n_obs <- 100
  x <- rnorm(n_obs)
  alpha <- 0
  beta1 <- 1
  eta1 <- alpha + as.matrix(x) %*% (beta1)
  s <-  0.1
  h_x_beta <- S(eta1,s)
  y <- rbinom(n = n_obs,size = 1,prob = h_x_beta)  
  dat2 <- data.frame(y,x,n=rnorm(n_obs))
  dat2$y <- as.factor(dat2$y)
  
  n_out <- 0
  outliers <- data.frame(x=rnorm(n_out,2),n=rnorm(n_out,2),y=rep(0,n_out))
  dat3 <- rbind(dat2,outliers)
  
  mod_log <- glm(formula = y ~ x , data = dat3, family = binomial(link = "logit"))
  mod_stu <- glm(formula = y ~ x , data = dat3, family = binomial(link = Gosset(0.8)))
  
  A1 <- as.data.frame(-coef(mod_log)[1]/coef(mod_log)[2])
  vlines <- data.frame(xint = c(-coef(mod_log)[1]/coef(mod_log)[2],-coef(mod_stu)[1]/coef(mod_stu)[2]),grp = c("logit", "st"))
  
  # plot1_po_pn <- 
  str(dat3)
  
  tikz(file = "img/Noisy_sou_0.tex", standAlone=F,width = 5.5, height = 2.8)
  # Correr todo para que funcione
  ggplot()+
    geom_point(dat3,mapping=aes(x=x,y=0,color=factor(y)),show.legend = F)+
    geom_vline(aes(xintercept = c(
      -coef(mod_log)[1]/coef(mod_log)[2],
      -coef(mod_stu)[1]/coef(mod_stu)[2]), 
      colour = c("#E69F00","#00AFBB")))+
    scale_colour_manual(values = c("#E69F00","#00AFBB","#B4464B", "#B4AF46"),
      labels = c("Logit","Student(0.8)", "",""))+
    labs(x = "x1", y = "", color = "")+
    theme(strip.background = element_rect(fill="yellow"),
      panel.background = element_rect(fill = "yellow", colour = "grey50"),
      plot.background = element_rect(fill = "yellow", colour = "yellow"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(fill = NA),
      legend.position = "top",
      legend.background = element_rect(fill = "yellow"))+
    # xlim(-3,3)+ylim(0,1)+
    scale_y_continuous(labels = NULL, breaks = NULL)
  endoffile <- dev.off()
  
  plot1_po_pn
  
  mod_log1 <- glm(formula = y ~ x + n, data = dat3, family = binomial(link = "logit"))
  mod_stu1 <- glm(formula = y ~ x + n, data = dat3, family = binomial(link = Gosset(0.8)))
  
  
  funcion1 <- function(x) (-coef(mod_log1)[1]/coef(mod_log1)[3])-((coef(mod_log1)[2]/coef(mod_log1)[3])*x)
  funcion1(-0.5)
  funcion1(0.5)
  
  funcion2 <- function(x) (-coef(mod_stu1)[1]/coef(mod_stu1)[3])-((coef(mod_stu1)[2]/coef(mod_stu1)[3])*x)
  funcion2(-0.5)
  funcion2(0.5)
  
  tikz(file = "img/Noisy_sou_2.tex", standAlone=F,width = 5.5, height = 2.8)
  ggplot()+
    geom_point(dat3,mapping=aes(x=x,y=as.numeric(as.character(n)),color=factor(y)),show.legend = F)+
    geom_segment(aes(y = c(
      funcion1(-2),
      funcion2(-2)),
      x = c(-2,-2),
      yend = c(
        funcion1(2),
        funcion2(2)
      ),
      xend = c(2,2),
      colour = c("#00AFBB","#E69F00")))+
    coord_cartesian(xlim = c(-3, 3), ylim = c(-2.5, 2.5))+
    scale_colour_manual(values = c("#E69F00","#00AFBB","#B4464B", "#B4AF46"),
      labels = c("Logit","Student(0.8)", "",""))+
    labs(x = "x1", y = "x2", color = "")+
    theme(strip.background = element_rect(fill="yellow"),
      panel.background = element_rect(fill = "yellow", colour = "grey50"),
      plot.background = element_rect(fill = "yellow", colour = "yellow"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.background = element_rect(fill = "yellow"),
      # legend.position = c(0.8,0.28),
      legend.key = element_rect(fill = NA)
      )
  endoffile <- dev.off()
  
  plot1_po_n <- ggplot(dat3,aes(x=x,y=n,color=y))+
    geom_point(show.legend = F)+
    scale_color_manual("", values = c("#B4464B", "#B4AF46"))+
    # stat_function(data = data.frame(x = seq(-1,1,0.001)), aes(x), inherit.aes = F,
    stat_function(fun = function(x) (-coef(mod_log1)[1]/coef(mod_log1)[3])-((coef(mod_log1)[2]/coef(mod_log1)[3])*x),
      aes(linetype = "dashed"), color = "black")+
    stat_function(fun = function(x) (-coef(mod_stu1)[1]/coef(mod_stu1)[3])-((coef(mod_stu1)[2]/coef(mod_stu1)[3])*x),
      aes(linetype = "solid"), color = "black")+
    ylim(-3,3)+
    xlim(-3,3)+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
      legend.key = element_rect(fill = "white"),
      # legend.position = c(2.1,-2))+
      legend.position = "top")+
    scale_linetype_manual(values = c("solid", "dashed"),
      labels = c("Student(0.8)", "Logit"))+
    labs(x = "x1", y = "x2", linetype  = "")
  
  plot1_po_n
  
  fun_log1 <- function(x) (-coef(mod_log1)[1]/coef(mod_log1)[3])-((coef(mod_log1)[2]/coef(mod_log1)[3])*x)
  fun_log_in1 <- function(x) (-x+(-coef(mod_log1)[1]/coef(mod_log1)[3]))/((coef(mod_log1)[2]/coef(mod_log1)[3]))
  fun_stu1 <- function(x) (-coef(mod_stu1)[1]/coef(mod_stu1)[3])-((coef(mod_stu1)[2]/coef(mod_stu1)[3])*x)
  fun_stu_in1 <- function(x) (-x+(-coef(mod_stu1)[1]/coef(mod_stu1)[3]))/((coef(mod_stu1)[2]/coef(mod_stu1)[3]))
  
  plot1_po_n <- ggplot(dat3,aes(x=x,y=n,color=y))+
    geom_point()+
    scale_color_manual("", values = c("#B4464B", "#B4AF46"))+
    geom_segment(color="black",linetype="solid",
      mapping = aes(x=fun_log_in1(-3),xend=fun_log_in1(3),y=fun_log1(fun_log_in1(-3)),yend=fun_log1(fun_log_in1(3))))+
    geom_segment(color="black",linetype="dashed",
      mapping = aes(x=fun_stu_in1(-3),xend=fun_stu_in1(3),y=fun_stu1(fun_stu_in1(-3)),yend=fun_stu1(fun_stu_in1(3))))+
    ylim(-3,3)+
    xlim(-3,3)+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
      legend.key = element_rect(fill = "white"),
      legend.position = "none")+
    ylab("$x_2$")+xlab("$x_1$")
  plot1_po_n
  
  n_out <- 5
  set.seed(seedd)
  outliers <- data.frame(x=rnorm(n_out,2,0.5),n=rnorm(n_out,0,0.5),y=rep(0,n_out))
  dat3 <- rbind(dat2,outliers)
  
  mod_log_n11 <- glm(formula = y ~ x + n, data = dat3, family = binomial(link = "logit"))
  mod_stu_n11 <- glm(formula = y ~ x + n , data = dat3, family = binomial(link = Gosset(0.8)))
  
  fun_log2 <- function(x) (-coef(mod_log_n11)[1]/coef(mod_log_n11)[3])-((coef(mod_log_n11)[2]/coef(mod_log_n11)[3])*x)
  fun_log_in2 <- function(x) (-x+(-coef(mod_log_n11)[1]/coef(mod_log_n11)[3]))/((coef(mod_log_n11)[2]/coef(mod_log_n11)[3]))
  fun_stu2 <- function(x) (-coef(mod_stu_n11)[1]/coef(mod_stu_n11)[3])-((coef(mod_stu_n11)[2]/coef(mod_stu_n11)[3])*x)
  fun_stu_in2 <- function(x) (-x+(-coef(mod_stu_n11)[1]/coef(mod_stu_n11)[3]))/((coef(mod_stu_n11)[2]/coef(mod_stu_n11)[3]))
  
  tikz(file = "img/Noisy_sou_3.tex", standAlone=F,width = 5.5, height = 2.8)
  ggplot()+
    geom_point(dat3,mapping=aes(x=x,y=as.numeric(as.character(n)),color=factor(y)),show.legend = F)+
    geom_segment(aes(y = c(
      fun_log2(fun_log_in2(-3)),
      fun_stu2(fun_stu_in2(-3))),
      x = c(fun_log_in2(-3),fun_stu_in2(-3)),
      yend = c(
        fun_log2(fun_log_in2(3)),
        fun_stu2(fun_stu_in2(3))
      ),
      xend = c(fun_log_in2(3),fun_stu_in2(3)),
      colour = c("#00AFBB","#E69F00")))+
    coord_cartesian(xlim = c(-3, 3), ylim = c(-2.5, 2.5))+
    scale_colour_manual(values = c("#E69F00","#00AFBB","#B4464B", "#B4AF46"),
      labels = c("Logit","Student(0.8)", "",""))+
    labs(x = "x1", y = "x2", color = "")+
    theme(strip.background = element_rect(fill="yellow"),
      panel.background = element_rect(fill = "yellow", colour = "grey50"),
      plot.background = element_rect(fill = "yellow", colour = "yellow"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.background = element_rect(fill = "yellow"),
      legend.key = element_rect(fill = NA))
  endoffile <- dev.off()
  
  plot1_ol_n <- ggplot(dat3,aes(x=x,y=n,color=y))+
    geom_point()+
    scale_color_manual("", values = c("#B4464B", "#B4AF46"))+
    geom_segment(color="black",linetype="solid",
      mapping = aes(x=fun_log_in2(-3),xend=fun_log_in2(3),y=fun_log2(fun_log_in2(-3)),yend=fun_log2(fun_log_in2(3))))+
    geom_segment(color="black",linetype="dashed",
      mapping = aes(x=fun_stu_in2(-3),xend=fun_stu_in2(3),y=fun_stu2(fun_stu_in2(-3)),yend=fun_stu2(fun_stu_in2(3))))+
    ylim(-3,3)+
    xlim(-3,3)+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
      legend.key = element_rect(fill = "white"),
      legend.position = "none")+
    ylab("$x_2$")+xlab("$x_1$")
  plot1_ol_n
  
  n_out <- 5
  set.seed(seedd)
  outliers <- data.frame(x=rnorm(n_out,2,0.5),n=rnorm(n_out,2,0.5),y=rep(0,n_out))
  dat3 <- rbind(dat2,outliers)
  
  mod_log_n1 <- glm(formula = y ~ x + n, data = dat3, family = binomial(link = "logit"))
  mod_stu_n1 <- glm(formula = y ~ x + n , data = dat3, family = binomial(link = Gosset(0.8)))
  
  fun_log3 <- function(x) (-coef(mod_log_n1)[1]/coef(mod_log_n1)[3])-((coef(mod_log_n1)[2]/coef(mod_log_n1)[3])*x)
  fun_log_in3 <- function(x) (-x+(-coef(mod_log_n1)[1]/coef(mod_log_n1)[3]))/((coef(mod_log_n1)[2]/coef(mod_log_n1)[3]))
  fun_stu3 <- function(x) (-coef(mod_stu_n1)[1]/coef(mod_stu_n1)[3])-((coef(mod_stu_n1)[2]/coef(mod_stu_n1)[3])*x)
  fun_stu_in3 <- function(x) (-x+(-coef(mod_stu_n1)[1]/coef(mod_stu_n1)[3]))/((coef(mod_stu_n1)[2]/coef(mod_stu_n1)[3]))
  
  tikz(file = "img/Noisy_sou_4.tex", standAlone=F,width = 5.5, height = 2.8)
  ggplot()+
    geom_point(dat3,mapping=aes(x=x,y=as.numeric(as.character(n)),color=factor(y)),show.legend = F)+
    geom_segment(aes(y = c(
      fun_log3(fun_log_in3(-3)),
      fun_stu3(fun_stu_in3(-3))),
      x = c(fun_log_in3(-3),fun_stu_in3(-3)),
      yend = c(
        fun_log3(fun_log_in3(3)),
        fun_stu3(fun_stu_in3(3))
      ),
      xend = c(fun_log_in3(3),fun_stu_in3(3)),
      colour = c("#00AFBB","#E69F00")))+
    coord_cartesian(xlim = c(-3, 3), ylim = c(-2, 2))+
    scale_colour_manual(values = c("#E69F00","#00AFBB","#B4464B", "#B4AF46"),
      labels = c("Logit","Student(0.8)", "",""))+
    labs(x = "x1", y = "x2", color = "")+
    theme(strip.background = element_rect(fill="yellow"),
      panel.background = element_rect(fill = "yellow", colour = "grey50"),
      plot.background = element_rect(fill = "yellow", colour = "yellow"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.background = element_rect(fill = "yellow"),
      legend.key = element_rect(fill = NA))
  endoffile <- dev.off()
  
  plot1_oa_n <- ggplot(dat3,aes(x=x,y=n,color=y))+
    geom_point()+
    scale_color_manual("", values = c("#B4464B", "#B4AF46"))+
    geom_segment(color="black",linetype="solid",
      mapping = aes(x=fun_log_in3(-3),xend=fun_log_in3(3),y=fun_log3(fun_log_in3(-3)),yend=fun_log3(fun_log_in3(3))))+
    geom_segment(color="black",linetype="dashed",
      mapping = aes(x=fun_stu_in3(-3),xend=fun_stu_in3(3),y=fun_stu3(fun_stu_in3(-3)),yend=fun_stu3(fun_stu_in3(3))))+
    ylim(-3,3)+
    xlim(-3,3)+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
      legend.key = element_rect(fill = "white"),
      legend.position = "none")+
    ylab("$x_2$")+xlab("$x_1$")
  plot1_oa_n
  
  
  ggarrange(plot1_po_pn,plot1_po_n,plot1_ol_n,plot1_oa_n,ncol = 2,nrow = 2)
}

tikz(file = "Numerical Example Noisy Variables.tex", standAlone=F,width = 7, height = 5)
# endoffile <- dev.off() 

set.seed(17)
n_obs <- 100
x <- rnorm(n_obs)
alpha <- 0
beta1 <- 1
eta1 <- alpha + as.matrix(x) %*% (beta1)
s <-  0.1
h_x_beta <- S(eta1,s)
y <- rbinom(n = n_obs,size = 1,prob = h_x_beta)  
dat2 <- data.frame(y,x,n=rnorm(n_obs))
dat2$y <- as.factor(dat2$y)

n_out <- 2
outliers <- data.frame(x=rnorm(n_out,2.5),n=rnorm(n_out,0),y=rep(0,n_out))
dat3 <- rbind(dat2,outliers)
head(dat3)
dat3$out <- 0

dat3[dat3$x == sort(dat3$x, decreasing = T)[1:2],"out"] <- 1
summary(dat3)
dat3$out <- as.factor(dat3$out)

plot1_ol_n_sa <- ggplot(dat3,aes(x=x,y=n,color=y,shape = out))+
  geom_point()+
  scale_color_manual("", values = c("#B4464B", "#B4AF46"))+
  # geom_segment(color="black",linetype="solid",
  #   mapping = aes(x=fun_log_in2(-3),xend=fun_log_in2(3),y=fun_log2(fun_log_in2(-3)),yend=fun_log2(fun_log_in2(3))))+
  # geom_segment(color="black",linetype="dashed",
  #   mapping = aes(x=fun_stu_in2(-3),xend=fun_stu_in2(3),y=fun_stu2(fun_stu_in2(-3)),yend=fun_stu2(fun_stu_in2(3))))+
  # ylim(-3,3)+
  # xlim(-2.7,2.3)+
  scale_shape_manual("", labels = c("0", "1"), values = c(19, 8))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none")+
  ylab("$x_2$")+xlab("$x_1$")
plot1_ol_n_sa

tikz(file = "tex/State of the art/outliers_ilustration.tex", standAlone=F,width = 5, height = 3)
# endoffile <- dev.off() 


summary(mod_log)
summary(mod_stu)
logLik(mod_stu)
BIC(mod_stu)
summary(mod_log_n1)
summary(mod_stu_n1)
logLik(mod_stu_n1)
BIC(mod_stu_n1)