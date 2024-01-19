# Numerical illustration --------------------------------------------------

plot1_po_pn <- plot1_po_n <- plot1_o_pn <- plot1_o_n <- NULL
set.seed(30) # 4  ou 11

num_exp_8a <- function(n_obs,s,n_out,df_stu,seeed,out_y="meth2"){
  # n_obs <- 200
  set.seed(seeed)
  x <- rnorm(n_obs)
  alpha <- 0
  # beta1 <- 1
  # beta1 <- c(.8, -.6, rep(0, 2))
  beta1 <- c(.8, .4, -.4, .2, rep(0, 4))
  # beta1 <- c(1, rep(0, 1))
  Sigma <- matrix(0, nrow = 8, ncol = 8)
  diag(Sigma) <- 1
  Ex <- cbind(MASS::mvrnorm(n_obs, rep(0, 8), Sigma))
  
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
  
  top3 <- order(Ex[,7],decreasing = T) # mayor a menor
  top3 <- dat2[top3,]
  top3 <- top3[top3$y==1,]
  
  top4 <- order(Ex[,8],decreasing = T) # mayor a menor
  top4 <- dat2[top4,]
  top4 <- top4[top4$y==1,]
  
  dat2[rownames(top)[1:round(n_out_f/4,digits = 0)],"y"] <- 0
  dat2[rownames(top2)[1:round(n_out_f/4,digits = 0)],"y"] <- 0
  dat2[rownames(top3)[1:round(n_out_f/4,digits = 0)],"y"] <- 0
  dat2[rownames(top4)[1:round(n_out_f/4,digits = 0)],"y"] <- 0
  
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
A8coefs <- A8AICsteplog <- A8BICstepstu <- A8AICstepstu <-  A8BICsteplog <- NULL
# for (outlier in 0:20) {
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp_8a(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable)
    A8coefs <- rbind(A8coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A8AICsteplog <- rbind(A8AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A8AICstepstu <- rbind(A8AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A8BICsteplog <- rbind(A8BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A8BICstepstu <- rbind(A8BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
  }
  print(outlier)
}

A8BICsteplog <- as.data.frame(A8BICsteplog)
A8BICstepstu <- as.data.frame(A8BICstepstu)
colnames(A8BICsteplog) <- c("BIC","outlier")
colnames(A8BICstepstu) <- c("BIC","outlier")
A8BICsteplog$link <- "logit"
A8BICstepstu$link <- "Student"
A8BICstep <- rbind(A8BICsteplog,A8BICstepstu)

A8BIC_PLOT <- ggplot(data = A8BICstep, aes(x=factor(outlier),y=BIC,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.05")

A8BIC_PLOT

head(A8coefs)
dat5 <- A8coefs[,-ncol(A8coefs)]
head(dat5)
dat5[dat5==0] <- NA
A8coefs[,-ncol(A8coefs)] <- dat5
colnames(A8coefs)[ncol(A8coefs)] <- "out"
A8summ_coef <- A8coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A8summ_coef)

A8_SUM <- cbind(A8summ_coef[,1:2],
  exva = rowMeans(A8summ_coef[,4:7]),
  novar = rowMeans(A8summ_coef[,8:11]))

A8BIC_PLOT <- ggplot(data = A8_SUM, aes(x=out,color=link)) +
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

A8BIC_PLOT

# 01 ---------------------------------------------------------------------

nsim <- 100
A81coefs <- A81AICsteplog <- A81BICstepstu <- A81AICstepstu <-  A81BICsteplog <- NULL
# for (outlier in 0:20) {
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp_8a(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable)
    A81coefs <- rbind(A81coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A81AICsteplog <- rbind(A81AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A81AICstepstu <- rbind(A81AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A81BICsteplog <- rbind(A81BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A81BICstepstu <- rbind(A81BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
  }
  print(outlier)
}

A81BICsteplog <- as.data.frame(A81BICsteplog)
A81BICstepstu <- as.data.frame(A81BICstepstu)
colnames(A81BICsteplog) <- c("BIC","outlier")
colnames(A81BICstepstu) <- c("BIC","outlier")
A81BICsteplog$link <- "logit"
A81BICstepstu$link <- "Student"
A81BICstep <- rbind(A81BICsteplog,A81BICstepstu)

A81BIC_PLOT <- ggplot(data = A81BICstep, aes(x=factor(outlier),y=BIC,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.1")

A81BIC_PLOT

head(A81coefs)
dat5 <- A81coefs[,-ncol(A81coefs)]
head(dat5)
dat5[dat5==0] <- NA
A81coefs[,-ncol(A81coefs)] <- dat5
colnames(A81coefs)[ncol(A81coefs)] <- "out"
A81summ_coef <- A81coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A81summ_coef)

A81_SUM <- cbind(A81summ_coef[,1:2],
  exva = rowMeans(A81summ_coef[,4:7]),
  novar = rowMeans(A81summ_coef[,8:11]))

A81BIC_PLOT <- ggplot(data = A81_SUM, aes(x=out,color=link)) +
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

A81BIC_PLOT

# 03 ---------------------------------------------------------------------

nsim <- 100
A83coefs <- A83AICsteplog <- A83BICstepstu <- A83AICstepstu <-  A83BICsteplog <- NULL
# for (outlier in 0:20) {
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp_8a(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable)
    A83coefs <- rbind(A83coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A83AICsteplog <- rbind(A83AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A83AICstepstu <- rbind(A83AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A83BICsteplog <- rbind(A83BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A83BICstepstu <- rbind(A83BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
  }
  print(outlier)
}

A83BICsteplog <- as.data.frame(A83BICsteplog)
A83BICstepstu <- as.data.frame(A83BICstepstu)
colnames(A83BICsteplog) <- c("BIC","outlier")
colnames(A83BICstepstu) <- c("BIC","outlier")
A83BICsteplog$link <- "logit"
A83BICstepstu$link <- "Student"
A83BICstep <- rbind(A83BICsteplog,A83BICstepstu)

A83BIC_PLOT <- ggplot(data = A83BICstep, aes(x=factor(outlier),y=BIC,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.3")

A83BIC_PLOT

head(A83coefs)
dat5 <- A83coefs[,-ncol(A83coefs)]
head(dat5)
dat5[dat5==0] <- NA
A83coefs[,-ncol(A83coefs)] <- dat5
colnames(A83coefs)[ncol(A83coefs)] <- "out"
A83summ_coef <- A83coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A83summ_coef)

A83_SUM <- cbind(A83summ_coef[,1:2],
  exva = rowMeans(A83summ_coef[,4:7]),
  novar = rowMeans(A83summ_coef[,8:11]))

A83BIC_PLOT <- ggplot(data = A83_SUM, aes(x=out,color=link)) +
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

A83BIC_PLOT

# 08 ---------------------------------------------------------------------

nsim <- 100
A88coefs <- A88AICsteplog <- A88BICstepstu <- A88AICstepstu <-  A88BICsteplog <- NULL
# for (outlier in 0:20) {
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp_8a(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable)
    A88coefs <- rbind(A88coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A88AICsteplog <- rbind(A88AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A88AICstepstu <- rbind(A88AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A88BICsteplog <- rbind(A88BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A88BICstepstu <- rbind(A88BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
  }
  print(outlier)
}

A88BICsteplog <- as.data.frame(A88BICsteplog)
A88BICstepstu <- as.data.frame(A88BICstepstu)
colnames(A88BICsteplog) <- c("BIC","outlier")
colnames(A88BICstepstu) <- c("BIC","outlier")
A88BICsteplog$link <- "logit"
A88BICstepstu$link <- "Student"
A88BICstep <- rbind(A88BICsteplog,A88BICstepstu)

A88BIC_PLOT <- ggplot(data = A88BICstep, aes(x=factor(outlier),y=BIC,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.05")

A88BIC_PLOT

head(A88coefs)
dat5 <- A88coefs[,-ncol(A88coefs)]
head(dat5)
dat5[dat5==0] <- NA
A88coefs[,-ncol(A88coefs)] <- dat5
colnames(A88coefs)[ncol(A88coefs)] <- "out"
A88summ_coef <- A88coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A88summ_coef)

A88_SUM <- cbind(A88summ_coef[,1:2],
  exva = rowMeans(A88summ_coef[,4:7]),
  novar = rowMeans(A88summ_coef[,8:11]))

A88BIC_PLOT <- ggplot(data = A88_SUM, aes(x=out,color=link)) +
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

A88BIC_PLOT



# Merge final -------------------------------------------------------------

library(ggpubr)


ggarrange(A8BIC_PLOT, A81BIC_PLOT, A83BIC_PLOT, A88BIC_PLOT, ncol = 2, nrow = 2) 


