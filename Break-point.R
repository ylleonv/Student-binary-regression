Gosset <- function(nu) {
  # qqt <- function(p, nu)
  #   sign(p-0.5)*sqrt(qf(1-2*pmin(p,1-p), 1, nu))
  # linkfun <- function(mu) qqt(mu,nu)
  linkfun <- function(mu) qt(mu, nu)
  linkinv <- function(eta) {
    # thresh <- -qqt(.Machine$double.eps,nu)
    # eta <- pmin(thresh, pmax(eta, -thresh))
    pt(eta, nu)
  }
  mu.eta <- function(eta) {
    pmax(dt(eta, nu), .Machine$double.eps)
  }
  valideta <- function(eta) TRUE
  name <- "Gosset"
  structure(list(
    linkfun = linkfun, linkinv = linkinv,
    mu.eta = mu.eta, valideta = valideta, name = name
  ),
    class = "link-glm"
  )
}

library(catdata)
library(ggplot2)
library(tidyr)
library(latex2exp)
library(tikzDevice)
library(ggpubr)
library(dplyr)

data(vaso)
vaso
vaso$vaso[vaso$vaso==2] <- 0

ggplot(data = vaso, aes(x = vol, y = rate, shape = factor(vaso))) +
  geom_point(aes(color = factor(vaso)))

response_plot <- ggplot(data = vaso, aes(x = vol, y = rate)) +
  # geom_point(aes(color = factor(vaso)))+
  # geom_point()+
  geom_text(
    label=vaso$vaso,
    # nudge_x = 0.25, nudge_y = 0.25,
    check_overlap = T
  )+
  geom_abline(intercept = 0.77, slope = -1.1, linetype="dashed",colour = "gray")+
  geom_abline(slope = 1, intercept = 0, color = "skyblue")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  geom_point(aes(y=link,x=vol), colour = "red")+
  xlab("ln of Volume") + ylab("ln of Rate")
response_plot

# tikz(file = "scatter_vaso.tex", standAlone=F,width = 3.5, height = 3)
response_plot
# endoffile <- dev.off() 


mod_0 <- glm(factor(vaso)~vol+rate, data = vaso,
  family = binomial(link = "logit"))
# 
coef(mod_0)
pred = -2.87 + 5.17*(seq(-1,1.5,by=0.01))+4.56*(seq(-1,1.5,by=0.01))
link <- predict(mod_0, type = "link")

plot_ly(x=vaso$vol, y=vaso$rate, z=link, type="scatter3d", mode="markers", color=vaso$vaso)


library(plotly)
library(tidyr)
library(dplyr)

dat1 <- tbl_df( cbind(link = pred, response = vaso$vaso))
str(dat1)

ggplot(data = dat1)+
  geom_point(aes(x=link,y=response))


ggplot(data = vaso, aes(x = vol, y = rate, colour = factor(vaso))) +
  # geom_point(aes(color = factor(vaso)))+
  geom_point()+
  # geom_text(
  #   label=vaso$vaso,
  #   # rownames(vaso),
  #   # nudge_x = 0.25, nudge_y = 0.25,
  #   check_overlap = T
  # )+
  scale_colour_manual(values=c("#B4464B", "#B4AF46"))+
  # geom_point(aes(y=pred, x= vol))+
  # geom_abline(slope = 1, intercept = 0, color = "skyblue", linetype = "twodash",xlim = c(-1,1.5))+
  stat_function(fun = function(x) (-coef(mod_0)[1]/coef(mod_0)[3])-((coef(mod_0)[2]/coef(mod_0)[3])*x),
    linetype = "solid", color = "gray", xlim = c(-1,1.5))+
  stat_function(fun = function(x) (0+1*x),
    color = "gray", linetype = "twodash", xlim = c(-1,1.5))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = "none")+
  xlab("ln of Volume") + ylab("ln of Rate")+  xlim(-1,1.5)
# geom_abline(intercept = 0.7, slope = -1.2)

tikz(file = "Scatterplot vaso.tex", standAlone=F,width = 4, height = 3)

# endoffile <- dev.off() 
 
mod_0 <- glm(factor(vaso)~vol+rate, data = vaso,
  family = binomial(link = Gosset(1)))

coef(mod_0)

new_dat <- rbind(vaso, c(-5, -5, 1))

mod_1 <- glm(factor(vaso)~vol+rate,  
  data = new_dat,
  family = binomial(link = "logit"))


ggplot(data = vaso, aes(x = vol, y = rate)) +
  # geom_point(aes(color = factor(vaso)))+
  # geom_point()+
  geom_text(
    label=vaso$vaso,
    # nudge_x = 0.25, nudge_y = 0.25,
    check_overlap = T
  )+
  # geom_point(aes(y=pred, x= vol))+
  geom_abline(slope = 1, intercept = 0, color = "skyblue")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab("ln of Volume") + ylab("ln of Rate")+
  geom_abline(intercept = -coef(mod_1)[1], slope = -coef(mod_1)[2])

mod_1 <- glm(factor(vaso)~vol+rate,  
  data = new_dat,
  family = binomial(link = Gosset(0.8)))


ggplot(data = vaso, aes(x = vol, y = rate)) +
  # geom_point(aes(color = factor(vaso)))+
  # geom_point()+
  geom_text(
    label=vaso$vaso,
    # nudge_x = 0.25, nudge_y = 0.25,
    check_overlap = T
  )+
  # geom_point(aes(y=pred, x= vol))+
  geom_abline(slope = 1, intercept = 0, color = "skyblue")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab("ln of Volume") + ylab("ln of Rate")+
  geom_abline(intercept = -coef(mod_1)[1], slope = -coef(mod_1)[2])



pred <- ifelse(predict(mod_0, type = "response")<0.5,0,1)

sum(diag(table(vaso$vaso, pred)))/nrow(vaso)



# Shapes ------------------------------------------------------------------

mod_1 <- glm(factor(vaso)~vol,  
  data = vaso,
  family = binomial(link = Gosset(0.6)))

dat1 <- data.frame(vol = seq(-1,1.5,by=0.065))

mod_2 <- glm(factor(vaso)~vol,  
  data = vaso,
  family = binomial(link = "logit"))

ggplot(data = vaso) +
  # geom_point(aes(color = factor(vaso)))+
  # geom_point(y=predict(mod_1, newdata = dat1, type = "response"), x = dat1$vol)+
  geom_path(y=as.vector(predict(mod_1, newdata = dat1, type = "response")),
    x = as.vector(dat1$vol),color = "skyblue")+
  geom_path(y=as.vector(predict(mod_2, newdata = dat1, type = "response")),
    x = as.vector(dat1$vol), color = "red")+
  geom_text(aes(x = vol, y = vaso),
    label=vaso$vaso,
    # nudge_x = 0.25, nudge_y = 0.25,
    check_overlap = T
  )+
  # geom_abline(slope = 1, intercept = 0, color = "skyblue")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab("ln of Volume") 

new_vaso2 <- rbind(vaso, c(5,0,1))

mod_1 <- glm(factor(vaso)~vol,  
  data = new_vaso2,
  family = binomial(link = Gosset(0.6)))

dat1 <- data.frame(vol = seq(-1,1.5,by=0.064))

mod_2 <- glm(factor(vaso)~vol,  
  data = new_vaso2,
  family = binomial(link = "logit"))

ggplot(data = new_vaso2) +
  # geom_point(aes(color = factor(vaso)))+
  # geom_point(y=predict(mod_1, newdata = dat1, type = "response"), x = dat1$vol)+
  geom_path(y=as.vector(predict(mod_1, newdata = dat1, type = "response")),
    x = as.vector(dat1$vol),color = "skyblue")+
  geom_path(y=as.vector(predict(mod_2, newdata = dat1, type = "response")),
    x = as.vector(dat1$vol), color = "red")+
  geom_text(aes(x = vol, y = vaso),
    label=new_vaso2$vaso,
    # nudge_x = 0.25, nudge_y = 0.25,
    check_overlap = T
  )+
  # geom_abline(slope = 1, intercept = 0, color = "skyblue")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab("ln of Volume") 


# Outlier (s,s,1)

range_s <- seq(-12,5,by = 0.5)

coef_ma <- matrix(nrow = length(range_s), ncol = 3)
accuracy <- matrix(nrow = length(range_s), ncol = 1)
ind <- 0


for (ind in 1:length(range_s)) {
  
  # ind <- 1
  new_dat <- rbind(vaso, c(range_s[ind], range_s[ind], 1))
  
  mod_1 <- glm(factor(vaso)~vol+rate,  
    data = new_dat,
    family = binomial(link = "logit"))
  
  coef_ma[ind,] <- coef(mod_1)
  
  pred <- ifelse(predict(mod_1, type = "response")<0.5,0,1)
  
  accuracy[ind] <- sum(diag(table(new_dat$vaso, pred)))/nrow(new_dat)
}

plot(accuracy)
plot(coef_ma[,2])

acc_dat <- data.frame(range_s, accuracy)

acc_plot <- ggplot(data = acc_dat) +
  geom_line(aes(x = range_s, y = accuracy))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab("s") + ylab("percent of correctly classified")
acc_plot

coef_ma <- cbind(coef_ma,range_s)
colnames(coef_ma) <- c("Intercept", "Beta 1", "Beta 2", "s")
# colnames(coef_ma) <- c("Beta 1", "Beta 2", "s")
coef_ma <- coef_ma[,-1]
coef_ma <- as.data.frame(coef_ma)

coef_to_plot <- gather(coef_ma, Beta, value, "Beta 1":"Beta 2")

betas_plot <- ggplot(data = coef_to_plot) +
  geom_line(aes(x = s,y=value,color=Beta))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.position = c(0.85, 0.2),)+
  # xlab("s") + ylab("% of correctly classified")+
  labs(title="", x="s", y="$X(s)$") +
  scale_color_manual("", labels = c("$X_1$", "$X_2$"), values = c("#4E84C4", "#52854C")) 

betas_plot

# tikz(file = "orginal_data_response.tex", standAlone=F,width = 3.5, height = 3)
response_plot

# tikz(file = "orginal_data_info.tex", standAlone=F,width = 6.5, height = 3)

# endoffile <- dev.off() 

ggarrange(betas_plot, acc_plot, ncol = 2)
# Student

# range_s <- seq(-12,15,by = 0.5)

coef_ma_s <- matrix(nrow = length(range_s), ncol = 3)
accuracy_s <- matrix(nrow = length(range_s), ncol = 1)
ind <- 0

for (ind in 1:length(range_s)) {
  
  new_dat <- rbind(vaso, c(range_s[ind], range_s[ind], 1))
  
  mod_1 <- glm(factor(vaso)~vol+rate,  
    data = new_dat,
    family = binomial(link = Gosset(0.63)))
  
  coef_ma_s[ind,] <- coef(mod_1)
  
  pred <- ifelse(predict(mod_1, type = "response")<0.5,0,1)
  
  accuracy_s[ind] <- sum(diag(table(new_dat$vaso, pred)))/nrow(new_dat)
}

coef_ma <- cbind(coef_ma_s,range_s)
colnames(coef_ma) <- c("Intercept", "Beta 1", "Beta 2", "s")
coef_ma <- coef_ma[,-1]
coef_ma <- as.data.frame(coef_ma)

coef_to_plot_s <- gather(coef_ma, Beta, value, "Beta 1":"Beta 2")
coef_to_plot$link <- "logistic"
coef_to_plot_s$link <- "Student(0.6)"
coef_tot <- rbind(coef_to_plot, coef_to_plot_s)

betas_plot <- ggplot(data = coef_tot) +
  geom_line(aes(x = s,y=value,color=Beta,linetype = link))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #       legend.position = c(0.85, 0.5),)+
  # xlab("s") + ylab("% of correctly classified")+
  labs(title="", x="s", y="$X(s)$") +
  scale_color_manual("", labels = c("$X_1$", "$X_2$"), values = c("#4E84C4", "#52854C")) 

betas_plot

# tikz(file = "betas_link.tex", standAlone=F, width = 6, height = 3)

# endoffile <- dev.off() 


# Normalize betas student -------------------------------------------------
library(GLMcat)

new_dat <- rbind(vaso, c(range_s[ind], range_s[ind], 1))

mod_1 <- glm(factor(vaso)~vol+rate,  
  data = new_dat,
  family = binomial(link = Gosset(0.6)))
coef(mod_1)

new_dat$vaso <- as.factor(new_dat$vaso)

mod_2 <- glmcat(vaso~vol+rate,  
  data = new_dat,
  cdf = list("student",0.6), 
  normalization = 0.95)

mod_2$normalization_s0* mod_2$coefficients

coef_ma_s1 <- matrix(nrow = length(range_s), ncol = 3)
accuracy_s1 <- matrix(nrow = length(range_s), ncol = 1)
ind <- 0


for (ind in 1:length(range_s)) {
  new_dat <- rbind(vaso, c(range_s[ind], range_s[ind], 1))
  new_dat$vaso <- as.factor(new_dat$vaso)
  mod_1 <- glmcat(vaso~vol+rate,  
    data = new_dat,
    cdf = list("student",0.8), 
    normalization = 0.95)
  
  coef_ma_s1[ind,] <- -mod_1$normalization_s0* mod_1$coefficients
  
  pred <- ifelse(predict(mod_1, type = "prob")<0.5,0,1)
  
  accuracy_s1[ind] <- sum(diag(table(new_dat$vaso, pred[,1])))/nrow(new_dat)
}

coef_ma <- cbind(coef_ma_s1,range_s)
colnames(coef_ma) <- c("Intercept", "Beta 1", "Beta 2", "s")
coef_ma <- coef_ma[,-1]
coef_ma <- as.data.frame(coef_ma)

coef_to_plot_s <- gather(coef_ma, Beta, value, "Beta 1":"Beta 2")
coef_to_plot$link <- "logistic"
coef_to_plot_s$link <- "Student(0.6)"
coef_tot <- rbind(coef_to_plot, coef_to_plot_s)


coef_ma_s[ind,] <- coef(mod_1)

mod_2 <- glmcat(vaso~vol+rate,  
  data = new_dat,
  # cdf = list("student",0.6), 
  normalization = 0.6)

mod_2$coefficients
mod_2$normalization_s0* mod_2$coefficients

head(coef_tot)
coef_tot$State <- paste0(coef_tot$Beta," ", coef_tot$link)
coef_tot$State <- as.factor(coef_tot$State)
summary(coef_tot)

coef_tot$State<- ordered(coef_tot$State,
  levels =c("Beta 1 Student(0.6)", "Beta 2 Student(0.6)", "Beta 1 logistic", "Beta 2 logistic"))
  
betas <- ggplot(data = coef_tot) +
  geom_line(aes(x = s,y=value,color=State, linetype = State))+
  scale_linetype_manual(values = c("dashed", "dashed","solid","solid")) +
  scale_color_manual(values = c("#4E84C4", "#52854C","#4E84C4", "#52854C")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.position = c(0.005, 0.75),
    legend.title = element_blank(),
    # legend.position="top",
    legend.justification = "left",
    # legend.key = element_rect(fill = "white")
    legend.key = element_blank(),
    # legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0.6, "pt"),
    legend.spacing.y = unit(0.6, "pt")
  )
betas
# tikz(file = "betas_link_norm.tex", standAlone=F, width = 6, height = 3)

endoffile <- dev.off()


acc_dat_s <- data.frame(range_s, accuracy = accuracy_s)
acc_dat$Link <- "logistic"
acc_dat_s$Link <- "Student(0.6)"
acc_tot <- rbind(acc_dat, acc_dat_s)

acc_plot <- ggplot(data = acc_tot) +
  geom_line(aes(x = range_s, y = accuracy, linetype = Link))+
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.position = c(0.55, 0.25),
    # legend.position="top", 
    legend.justification = "left",
    legend.key = element_rect(fill = "white")
    # legend.margin = margin(0, 0, 0, 0),
    # legend.spacing.x = unit(1, "pt")
  )+
  guides(linetype = guide_legend(title = NULL)) +
  xlab("(a) s") + ylab("Accuracy")
acc_plot

library(ggpubr)
ggarrange(betas,acc_plot,ncol = 2, nrow = 1, labels = c("",""))


# tikz(file = "betas_acc.tex", standAlone=F, width = 7.5, height = 3)
# endoffile <- dev.off()



# Separation issue


mod_0 <- glm(factor(vaso)~vol+rate, data = vaso,
  family = binomial(link = "logit"))

coef(mod_0)
summary(mod_0)[,4]

mod_1 <- glm(factor(vaso)~vol+rate, data = vaso[-c(4,18),],
  family = binomial(link = "logit"))

coef(mod_1)
summary(mod_1)

mod_2 <- glm(factor(vaso)~vol+rate, data = vaso[-c(4,18),],
  family = binomial(link = Gosset(0.8)))

coef(mod_2)
summary(mod_2)


pred <- ifelse(predict(mod_0, type = "response")<0.5,0,1)

# Stability of the p value

vaso1 <- vaso[-c(4,18),]

ggplot(data = vaso1, aes(x = vol, y = rate)) +
  # geom_point(aes(color = factor(vaso)))+
  # geom_point()+
  geom_text(
    label=vaso1$vaso,
    # nudge_x = 0.25, nudge_y = 0.25,
    check_overlap = T
  )+
  geom_abline(slope = 1, intercept = 0, color = "skyblue")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab("ln of Volume") + ylab("ln of Rate")

range_s <- seq(-12,5,by = 0.5)

pval_ma <- pval_ma_s <-  matrix(nrow = length(range_s), ncol = 3)
accuracy <- matrix(nrow = length(range_s), ncol = 1)
ind <- 0


for (ind in 1:length(range_s)) {
  
  new_dat <- rbind(vaso1, c(range_s[ind], range_s[ind], 1))
  
  mod_1 <- glm(factor(vaso)~vol+rate,  
    data = new_dat,
    family = binomial(link = "logit"))
  dat1 <- summary(mod_1)
  
  dat1$coefficients[,4]
  
  pval_ma[ind,] <- dat1$coefficients[,4]
  
  # pred <- ifelse(predict(mod_1, type = "response")<0.5,0,1)
  # 
  # accuracy[ind] <- sum(diag(table(new_dat$vaso, pred)))/nrow(new_dat)
}



for (ind in 1:length(range_s)) {
  
  new_dat <- rbind(vaso1, c(range_s[ind], range_s[ind], 1))
  
  mod_1 <- glm(factor(vaso)~vol+rate,  
    data = new_dat,
    family = binomial(link = Gosset(0.8)))
  dat1 <- summary(mod_1)
  
  dat1$coefficients[,4]
  
  pval_ma_s[ind,] <- dat1$coefficients[,4]
  
  # pred <- ifelse(predict(mod_1, type = "response")<0.5,0,1)
  # 
  # accuracy[ind] <- sum(diag(table(new_dat$vaso, pred)))/nrow(new_dat)
}

intercept_p <- data.frame(intercept = pval_ma[,1], range_s)
intercept_p_s <- data.frame(intercept = pval_ma_s[,1], range_s)

intercept_p$Link <- "logistic"
intercept_p_s$Link <- "Student(0.8)"

p_val_int <- rbind(intercept_p,intercept_p_s)

intercept <- ggplot(data = p_val_b1)+
  geom_line(aes(x=range_s,y=intercept,linetype = Link))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.position = "none") +
  xlab("s") + ylab("p-value") +
  geom_hline(yintercept=0.05, linetype = "twodash", color = "skyblue")
intercept

# Change pval_ma_s[,1] by pval_ma_s[,2] above

intercept_p <- data.frame(intercept = pval_ma[,2], range_s)
intercept_p_s <- data.frame(intercept = pval_ma_s[,2], range_s)

intercept_p$Link <- "logistic"
intercept_p_s$Link <- "Student(0.8)"

p_val_b1 <- rbind(intercept_p,intercept_p_s)
beta1 <- ggplot(data = p_val_b1)+
  geom_line(aes(x=range_s,y=intercept,linetype = Link))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.position = "none") +
  xlab("s") + ylab("p-value")

intercept_p <- data.frame(intercept = pval_ma[,3], range_s)
intercept_p_s <- data.frame(intercept = pval_ma_s[,3], range_s)

intercept_p$Link <- "logistic"
intercept_p_s$Link <- "Student(0.8)"

p_val_b2 <- rbind(intercept_p,intercept_p_s)
beta2 <- ggplot(data = p_val_b1)+
  geom_line(aes(x=range_s,y=intercept,linetype = Link))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.position = "none") +
  xlab("s") + ylab("p-value")

p_val_int$Beta <- "0Intercept"
p_val_b1$Beta <- "beta1"
p_val_b2$Beta <- "beta2"

p_val_tot <- rbind(p_val_int,p_val_b1,p_val_b2)

ggplot(data = p_val_tot)+
  geom_line(aes(x=range_s,y=intercept,linetype = Link))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.position = "none") +
  xlab("s") + ylab("p-value") +
  geom_hline(yintercept=0.05, linetype = "twodash", color = "skyblue")+
  facet_grid(cols = vars(Beta))



# tikz(file = "p_stability.tex", standAlone=F,width = 6.5, height = 3)

# endoffile <- dev.off() 

