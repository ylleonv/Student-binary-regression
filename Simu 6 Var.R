# Numerical illustration --------------------------------------------------

plot1_po_pn <- plot1_po_n <- plot1_o_pn <- plot1_o_n <- NULL
set.seed(30) # 4  ou 11

num_exp_6a <- function(n_obs,s,n_out,df_stu,seeed,out_y="meth2"){
  # n_obs <- 200
  set.seed(seeed)
  x <- rnorm(n_obs)
  alpha <- 0
  # beta1 <- 1
  # beta1 <- c(.8, -.6, rep(0, 2))
  beta1 <- c(.8, .4, -.4, .2, rep(0, 2))
  # beta1 <- c(1, rep(0, 1))
  Sigma <- matrix(0, nrow = 6, ncol = 6)
  diag(Sigma) <- 1
  Ex <- cbind(MASS::mvrnorm(n_obs, rep(0, 6), Sigma))
  
  eta1 <- alpha + as.matrix(Ex) %*% (beta1)
  # s <-  0.05
  h_x_beta <- S(eta1,s)
  y <- rbinom(n = n_obs,size = 1,prob = h_x_beta)  
  dat2 <- data.frame(y,Ex)
  dat2$y <- as.factor(dat2$y)
  
  # n_out <- 5
  n_out_f <- n_out
  # top <- order(h_x_beta,decreasing = T) # mayor a menor
  top <- order(Ex[,5],decreasing = T) # mayor a menor
  top <- dat2[top,]
  top <- top[top$y==1,]
  
  top2 <- order(Ex[,6],decreasing = T) # mayor a menor
  top2 <- dat2[top2,]
  top2 <- top2[top2$y==1,]
  
  dat2[rownames(top)[1:round(n_out_f/2,digits = 0)],"y"] <- 0
  dat2[rownames(top2)[1:round(n_out_f/2,digits = 0)],"y"] <- 0
  # dat2[rownames(top3)[1:round(n_out_f/4,digits = 0)],"y"] <- 0
  # dat2[rownames(top4)[1:round(n_out_f/4,digits = 0)],"y"] <- 0
  
  dat3 <- dat2
  
  
  mod_log_n1 <- glm(formula = y ~ . , data = dat3, family = binomial(link = "logit"))
  mod_stu_n1 <- glm(formula = y ~ . , data = dat3, family = binomial(link = Gosset(df_stu)))
  
  step_log <- step(mod_log_n1, direction = "backward", trace = F, k=log(nrow(dat3)))
  step_stu <- step(mod_stu_n1, direction = "backward", trace = F, k=log(nrow(dat3)))
  
  step_log_c <- step_log$coefficients
  step_stu_c <- step_stu$coefficients
  
  pv_log <- summary(mod_log_n1)$coefficients[,4]
  
  dat_sep <- bind_rows(
    as.data.frame(cbind(t(pv_log),link = "BB")),
    as.data.frame(cbind(t(step_log_c),link = "logit")),
    as.data.frame(cbind(t(step_stu_c),link = "stu")))
  dat_sep <- dat_sep[-1,]
  dat_sep[is.na(dat_sep)] <- 0
  dat_b <- dat_sep
  
  # plot1_o_n <- ggplot(dat3,aes(x=x,y=n,color=y))+
  #   geom_point()+
  #   scale_color_manual("", values = c("#E69F00", "#00AFBB"))+
  #   stat_function(fun = function(x) (-coef(mod_log_n1)[1]/coef(mod_log_n1)[3])-((coef(mod_log_n1)[2]/coef(mod_log_n1)[3])*x),
  #     linetype = "solid", color = "black")+
  #   stat_function(fun = function(x) (-coef(mod_stu_n1)[1]/coef(mod_stu_n1)[3])-((coef(mod_stu_n1)[2]/coef(mod_stu_n1)[3])*x),
  #     linetype = "dashed", color = "black")+
  #   ylim(-3,3)+
  #   theme(panel.background = element_rect(fill = "white", colour = "grey50"),
  #     legend.key = element_rect(fill = "white"),
  #     legend.position = "none")
  # 
  # print(ggarrange(plot1_po_pn,plot1_po_n,plot1_o_pn,plot1_o_n,ncol = 2,nrow = 2))
  
  return(list(
    # sum_log = summary(mod_log),
    # sum_stu = summary(mod_stu),
    sum_log_1 = summary(mod_log_n1),
    sum_stu_1 = summary(mod_stu_n1),
    # log_log = logLik(mod_log),
    log_log1 = logLik(mod_log_n1),
    # log_stu = logLik(mod_stu),
    log_stu1 = logLik(mod_stu_n1),
    # bic_log = BIC(mod_log),
    bic_log1 = BIC(mod_log_n1),
    # bic_stu = BIC(mod_stu),
    bic_stu1 = BIC(mod_stu_n1),
    
    # AIC_log = AIC(mod_log),
    AIC_log1 = AIC(mod_log_n1),
    # AIC_stu = AIC(mod_stu),
    AIC_stu1 = AIC(mod_stu_n1),
    AIC_step_log = AIC(step_log),
    AIC_step_stu = AIC(step_stu),
    BIC_step_log = BIC(step_log),
    BIC_step_stu = BIC(step_stu),
    coef_sel = dat_b))
}

# 0.05 --------------------------------------------------------------------

nsim <- 100
A6coefs <- A6AICsteplog <- A6BICstepstu <- A6AICstepstu <-  A6BICsteplog <- NULL
# for (outlier in 0:20) {
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp_6a(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable)
    A6coefs <- rbind(A6coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A6AICsteplog <- rbind(A6AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A6AICstepstu <- rbind(A6AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A6BICsteplog <- rbind(A6BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A6BICstepstu <- rbind(A6BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
  }
  print(outlier)
}

A6BICsteplog <- as.data.frame(A6BICsteplog)
A6BICstepstu <- as.data.frame(A6BICstepstu)
colnames(A6BICsteplog) <- c("BIC","outlier")
colnames(A6BICstepstu) <- c("BIC","outlier")
A6BICsteplog$link <- "logit"
A6BICstepstu$link <- "Student"
A6BICstep <- rbind(A6BICsteplog,A6BICstepstu)

A6BIC_PLOT <- ggplot(data = A6BICstep, aes(x=factor(outlier),y=BIC,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.05")

A6BIC_PLOT

head(A6coefs)
dat5 <- A6coefs[,-ncol(A6coefs)]
head(dat5)
dat5[dat5==0] <- NA
A6coefs[,-ncol(A6coefs)] <- dat5
colnames(A6coefs)[ncol(A6coefs)] <- "out"
A6summ_coef <- A6coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A6summ_coef)

A6_SUM <- cbind(A6summ_coef[,1:2],
  exva = rowMeans(A6summ_coef[,4:7]),
  novar = rowMeans(A6summ_coef[,8:9]))

A6BIC_PLOT <- ggplot(data = A6_SUM, aes(x=out,color=link)) +
  geom_line(aes(y=novar, linetype = "dashed")) +
  geom_line(aes(y=exva,linetype = "solid")) +
  geom_point(aes(y=novar))+
  geom_point(aes(y=exva))+
  xlab("") + ylab("") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.05")

A6BIC_PLOT

# 01 ---------------------------------------------------------------------

nsim <- 100
A61coefs <- A61AICsteplog <- A61BICstepstu <- A61AICstepstu <-  A61BICsteplog <- NULL
# for (outlier in 0:20) {
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp_6a(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable)
    A61coefs <- rbind(A61coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A61AICsteplog <- rbind(A61AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A61AICstepstu <- rbind(A61AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A61BICsteplog <- rbind(A61BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A61BICstepstu <- rbind(A61BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
  }
  print(outlier)
}

A61BICsteplog <- as.data.frame(A61BICsteplog)
A61BICstepstu <- as.data.frame(A61BICstepstu)
colnames(A61BICsteplog) <- c("BIC","outlier")
colnames(A61BICstepstu) <- c("BIC","outlier")
A61BICsteplog$link <- "logit"
A61BICstepstu$link <- "Student"
A61BICstep <- rbind(A61BICsteplog,A61BICstepstu)

A61BIC_PLOT <- ggplot(data = A61BICstep, aes(x=factor(outlier),y=BIC,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.1")

A61BIC_PLOT

head(A61coefs)
dat5 <- A61coefs[,-ncol(A61coefs)]
head(dat5)
dat5[dat5==0] <- NA
A61coefs[,-ncol(A61coefs)] <- dat5
colnames(A61coefs)[ncol(A61coefs)] <- "out"
A61summ_coef <- A61coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A61summ_coef)

A61_SUM <- cbind(A61summ_coef[,1:2],
  exva = rowMeans(A61summ_coef[,4:7]),
  novar = rowMeans(A61summ_coef[,8:9]))

A61BIC_PLOT <- ggplot(data = A61_SUM, aes(x=out,color=link)) +
  geom_line(aes(y=novar, linetype = "dashed")) +
  geom_line(aes(y=exva,linetype = "solid")) +
  geom_point(aes(y=novar))+
  geom_point(aes(y=exva))+
  xlab("") + ylab("") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.1")

A61BIC_PLOT

# 03 ---------------------------------------------------------------------

nsim <- 100
A63coefs <- A63AICsteplog <- A63BICstepstu <- A63AICstepstu <-  A63BICsteplog <- NULL
# for (outlier in 0:20) {
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp_6a(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable)
    A63coefs <- rbind(A63coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A63AICsteplog <- rbind(A63AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A63AICstepstu <- rbind(A63AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A63BICsteplog <- rbind(A63BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A63BICstepstu <- rbind(A63BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
  }
  print(outlier)
}

A63BICsteplog <- as.data.frame(A63BICsteplog)
A63BICstepstu <- as.data.frame(A63BICstepstu)
colnames(A63BICsteplog) <- c("BIC","outlier")
colnames(A63BICstepstu) <- c("BIC","outlier")
A63BICsteplog$link <- "logit"
A63BICstepstu$link <- "Student"
A63BICstep <- rbind(A63BICsteplog,A63BICstepstu)

A63BIC_PLOT <- ggplot(data = A63BICstep, aes(x=factor(outlier),y=BIC,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.3")

A63BIC_PLOT

head(A63coefs)
dat5 <- A63coefs[,-ncol(A63coefs)]
head(dat5)
dat5[dat5==0] <- NA
A63coefs[,-ncol(A63coefs)] <- dat5
colnames(A63coefs)[ncol(A63coefs)] <- "out"
A63summ_coef <- A63coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A63summ_coef)

A63_SUM <- cbind(A63summ_coef[,1:2],
  exva = rowMeans(A63summ_coef[,4:7]),
  novar = rowMeans(A63summ_coef[,8:9]))

A63BIC_PLOT <- ggplot(data = A63_SUM, aes(x=out,color=link)) +
  geom_line(aes(y=novar, linetype = "dashed")) +
  geom_line(aes(y=exva,linetype = "solid")) +
  geom_point(aes(y=novar))+
  geom_point(aes(y=exva))+
  xlab("") + ylab("") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.3")

A63BIC_PLOT

# 08 ---------------------------------------------------------------------

nsim <- 100
A68coefs <- A68AICsteplog <- A68BICstepstu <- A68AICstepstu <-  A68BICsteplog <- NULL
# for (outlier in 0:20) {
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp_6a(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable)
    A68coefs <- rbind(A68coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A68AICsteplog <- rbind(A68AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A68AICstepstu <- rbind(A68AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A68BICsteplog <- rbind(A68BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A68BICstepstu <- rbind(A68BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
  }
  print(outlier)
}

A68BICsteplog <- as.data.frame(A68BICsteplog)
A68BICstepstu <- as.data.frame(A68BICstepstu)
colnames(A68BICsteplog) <- c("BIC","outlier")
colnames(A68BICstepstu) <- c("BIC","outlier")
A68BICsteplog$link <- "logit"
A68BICstepstu$link <- "Student"
A68BICstep <- rbind(A68BICsteplog,A68BICstepstu)

A68BIC_PLOT <- ggplot(data = A68BICstep, aes(x=factor(outlier),y=BIC,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.05")

A68BIC_PLOT

head(A68coefs)
dat5 <- A68coefs[,-ncol(A68coefs)]
head(dat5)
dat5[dat5==0] <- NA
A68coefs[,-ncol(A68coefs)] <- dat5
colnames(A68coefs)[ncol(A68coefs)] <- "out"
A68summ_coef <- A68coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A68summ_coef)

A68_SUM <- cbind(A68summ_coef[,1:2],
  exva = rowMeans(A68summ_coef[,4:7]),
  novar = rowMeans(A68summ_coef[,8:9]))

A68BIC_PLOT <- ggplot(data = A68_SUM, aes(x=out,color=link)) +
  geom_line(aes(y=novar, linetype = "dashed")) +
  geom_line(aes(y=exva,linetype = "solid")) +
  geom_point(aes(y=novar))+
  geom_point(aes(y=exva))+
  xlab("") + ylab("") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.8")

A68BIC_PLOT



# Merge final -------------------------------------------------------------

library(ggpubr)

ggarrange(A6BIC_PLOT, A61BIC_PLOT, A63BIC_PLOT, A68BIC_PLOT, ncol = 2, nrow = 2) 


