# Numerical illustration --------------------------------------------------

plot1_po_pn <- plot1_po_n <- plot1_o_pn <- plot1_o_n <- NULL
set.seed(30) # 4  ou 11

num_exp_2a <- function(n_obs,s,n_out,df_stu,seeed,out_y="meth2"){
  # n_obs <- 200
  set.seed(seeed)
  x <- rnorm(n_obs)
  alpha <- 0
  # beta1 <- 1
  beta1 <- c(.8, -.6, rep(0, 2))
  # beta1 <- c(1, rep(0, 1))
  Sigma <- matrix(0, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  Ex <- cbind(MASS::mvrnorm(n_obs, rep(0, 4), Sigma))
  
  eta1 <- alpha + as.matrix(Ex) %*% (beta1)
  # s <-  0.05
  h_x_beta <- S(eta1,s)
  y <- rbinom(n = n_obs,size = 1,prob = h_x_beta)  
  dat2 <- data.frame(y,Ex)
  dat2$y <- as.factor(dat2$y)
  
  # n_out <- 5
  n_out_f <- n_out
  # top <- order(h_x_beta,decreasing = T) # mayor a menor
  top <- order(Ex[,3],decreasing = T) # mayor a menor
  top <- dat2[top,]
  top <- top[top$y==1,]
  
  top2 <- order(Ex[,4],decreasing = T) # mayor a menor
  top2 <- dat2[top2,]
  top2 <- top2[top2$y==1,]
  
  dat2[rownames(top)[1:round(n_out_f/2,digits = 0)],"y"] <- 0
  dat2[rownames(top2)[1:round(n_out_f/2,digits = 0)],"y"] <- 0
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
A2coefs <- A2AICsteplog <- A2BICstepstu <- A2AICstepstu <-  A2BICsteplog <- NULL
# for (outlier in 0:20) {
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp_2a(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable)
    A2coefs <- rbind(A2coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A2AICsteplog <- rbind(A2AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A2AICstepstu <- rbind(A2AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A2BICsteplog <- rbind(A2BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A2BICstepstu <- rbind(A2BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
  }
  print(outlier)
}

A2BICsteplog <- as.data.frame(A2BICsteplog)
A2BICstepstu <- as.data.frame(A2BICstepstu)
colnames(A2BICsteplog) <- c("BIC","outlier")
colnames(A2BICstepstu) <- c("BIC","outlier")
A2BICsteplog$link <- "logit"
A2BICstepstu$link <- "Student"
A2BICstep <- rbind(A2BICsteplog,A2BICstepstu)

A2BIC_PLOT <- ggplot(data = A2BICstep, aes(x=factor(outlier),y=BIC,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.05")

A2BIC_PLOT

head(A2coefs)
dat5 <- A2coefs[,-ncol(A2coefs)]
head(dat5)
dat5[dat5==0] <- NA
A2coefs[,-ncol(A2coefs)] <- dat5
colnames(A2coefs)[ncol(A2coefs)] <- "out"
A2summ_coef <- A2coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A2summ_coef)

A2_SUM <- cbind(A2summ_coef[,1:2],
  exva = rowMeans(A2summ_coef[,4:5]),
  novar = rowMeans(A2summ_coef[,6:7]))

A2BIC_PLOT <- ggplot(data = A2_SUM, aes(x=out,color=link)) +
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

A2BIC_PLOT

# 01 ---------------------------------------------------------------------

nsim <- 100
A21coefs <- A21AICsteplog <- A21BICstepstu <- A21AICstepstu <-  A21BICsteplog <- NULL
# for (outlier in 0:20) {
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp_2a(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable)
    A21coefs <- rbind(A21coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A21AICsteplog <- rbind(A21AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A21AICstepstu <- rbind(A21AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A21BICsteplog <- rbind(A21BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A21BICstepstu <- rbind(A21BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
  }
  print(outlier)
}

A21BICsteplog <- as.data.frame(A21BICsteplog)
A21BICstepstu <- as.data.frame(A21BICstepstu)
colnames(A21BICsteplog) <- c("BIC","outlier")
colnames(A21BICstepstu) <- c("BIC","outlier")
A21BICsteplog$link <- "logit"
A21BICstepstu$link <- "Student"
A21BICstep <- rbind(A21BICsteplog,A21BICstepstu)

A21BIC_PLOT <- ggplot(data = A21BICstep, aes(x=factor(outlier),y=BIC,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.1")

A21BIC_PLOT

head(A21coefs)
dat5 <- A21coefs[,-ncol(A21coefs)]
head(dat5)
dat5[dat5==0] <- NA
A21coefs[,-ncol(A21coefs)] <- dat5
colnames(A21coefs)[ncol(A21coefs)] <- "out"
A21summ_coef <- A21coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A21summ_coef)

A21_SUM <- cbind(A21summ_coef[,1:2],
  exva = rowMeans(A21summ_coef[,4:5]),
  novar = rowMeans(A21summ_coef[,6:7]))

A21BIC_PLOT <- ggplot(data = A21_SUM, aes(x=out,color=link)) +
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

A21BIC_PLOT

# 03 ---------------------------------------------------------------------

nsim <- 100
A23coefs <- A23AICsteplog <- A23BICstepstu <- A23AICstepstu <-  A23BICsteplog <- NULL
# for (outlier in 0:20) {
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp_2a(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable)
    A23coefs <- rbind(A23coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A23AICsteplog <- rbind(A23AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A23AICstepstu <- rbind(A23AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A23BICsteplog <- rbind(A23BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A23BICstepstu <- rbind(A23BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
  }
  print(outlier)
}

A23BICsteplog <- as.data.frame(A23BICsteplog)
A23BICstepstu <- as.data.frame(A23BICstepstu)
colnames(A23BICsteplog) <- c("BIC","outlier")
colnames(A23BICstepstu) <- c("BIC","outlier")
A23BICsteplog$link <- "logit"
A23BICstepstu$link <- "Student"
A23BICstep <- rbind(A23BICsteplog,A23BICstepstu)

A23BIC_PLOT <- ggplot(data = A23BICstep, aes(x=factor(outlier),y=BIC,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.3")

A23BIC_PLOT

head(A23coefs)
dat5 <- A23coefs[,-ncol(A23coefs)]
head(dat5)
dat5[dat5==0] <- NA
A23coefs[,-ncol(A23coefs)] <- dat5
colnames(A23coefs)[ncol(A23coefs)] <- "out"
A23summ_coef <- A23coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A23summ_coef)

A23_SUM <- cbind(A23summ_coef[,1:2],
  exva = rowMeans(A23summ_coef[,4:5]),
  novar = rowMeans(A23summ_coef[,6:7]))

A23BIC_PLOT <- ggplot(data = A23_SUM, aes(x=out,color=link)) +
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

A23BIC_PLOT

# 08 ---------------------------------------------------------------------

nsim <- 100
A28coefs <- A28AICsteplog <- A28BICstepstu <- A28AICstepstu <-  A28BICsteplog <- NULL
# for (outlier in 0:20) {
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp_2a(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable)
    A28coefs <- rbind(A28coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A28AICsteplog <- rbind(A28AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A28AICstepstu <- rbind(A28AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A28BICsteplog <- rbind(A28BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A28BICstepstu <- rbind(A28BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
  }
  print(outlier)
}

A28BICsteplog <- as.data.frame(A28BICsteplog)
A28BICstepstu <- as.data.frame(A28BICstepstu)
colnames(A28BICsteplog) <- c("BIC","outlier")
colnames(A28BICstepstu) <- c("BIC","outlier")
A28BICsteplog$link <- "logit"
A28BICstepstu$link <- "Student"
A28BICstep <- rbind(A28BICsteplog,A28BICstepstu)

A28BIC_PLOT <- ggplot(data = A28BICstep, aes(x=factor(outlier),y=BIC,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.05")

A28BIC_PLOT

head(A28coefs)
dat5 <- A28coefs[,-ncol(A28coefs)]
head(dat5)
dat5[dat5==0] <- NA
A28coefs[,-ncol(A28coefs)] <- dat5
colnames(A28coefs)[ncol(A28coefs)] <- "out"
A28summ_coef <- A28coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A28summ_coef)

A28_SUM <- cbind(A28summ_coef[,1:2],
  exva = rowMeans(A28summ_coef[,4:5]),
  novar = rowMeans(A28summ_coef[,6:7]))

A28BIC_PLOT <- ggplot(data = A28_SUM, aes(x=out,color=link)) +
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

A28BIC_PLOT



# Merge final -------------------------------------------------------------

library(ggpubr)


ggarrange(A2BIC_PLOT, A21BIC_PLOT, A23BIC_PLOT, A28BIC_PLOT, ncol = 2, nrow = 2) 


