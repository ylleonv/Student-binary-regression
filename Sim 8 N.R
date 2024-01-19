library(ggplot2)
# Numerical illustration --------------------------------------------------

plot1_po_pn <- plot1_po_n <- plot1_o_pn <- plot1_o_n <- NULL
set.seed(30) # 4  ou 11
g_deviance <- function(nu, data_in) {
  logLik(glm(y~.,family = binomial(Gosset(nu)),data = data_in, maxit = 100))
}

num_exp_8 <- function(n_obs,s,n_out,df_stu,seeed,out_y="meth2"){
  set.seed(seeed)
  alpha <- 0
  # beta1 <- c(1, rep(0, 1))
  Sigma <- matrix(0, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  beta1 <- c(.8, 0.4, -0.4, 0.2, rep(0, 4))
  
  Sigma <- matrix(0, nrow = 8, ncol = 8)
  diag(Sigma) <- 1
  Ex <- cbind(MASS::mvrnorm(n_obs, rep(0, 8), Sigma))
  
  eta1 <- alpha + as.matrix(Ex) %*% (beta1)
  h_x_beta <- S(eta1,s)
  y <- rbinom(n = n_obs,size = 1,prob = h_x_beta)  
  dat2 <- data.frame(y,Ex)
  dat2$y <- as.factor(dat2$y)
  dat3 <- dat2
  # df_stu <- 0.8
  # mod_log <- glm(formula = y ~ X1 , data = dat3, family = binomial(link = "logit"))
  # mod_stu <- glm(formula = y ~ X1 , data = dat3, family = binomial(link = Gosset(df_stu)))
  
  n_out_f <- n_out
  Sigma2 <- matrix(0, nrow = 8, ncol = 8)
  diag(Sigma2) <- 0.25
  if(n_out_f > 0){
    if(n_out_f==1){
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = 2, c(1,1,-1,0.5,2,2,2,2), Sigma2)),y=rep(0,n_out_f)) 
      dat3 <- rbind(dat2,outliers[1,])
    }else{
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = n_out_f, c(1,1,-1,0.5,2,2,2,2), Sigma2)),y=rep(0,n_out_f)) 
      dat3 <- rbind(dat2,outliers)
    }
  }
  
  nu_expl <- c(seq(0.3,1,0.1),2:8)
  log_nu <- NULL
  i <- 1
  for (nu in nu_expl) {
    log_nu[i] <- logLik(glm(formula = y ~ ., data = dat3, family = binomial(link = Gosset(nu))))
    i <- i+1
  }
  max_grille <- which.max(log_nu)
  max_grille <- nu_expl[max_grille]
  result_nu <- optimize(g_deviance, c(0.25, 1), maximum = T, dat3)
  df_stu_opt <- result_nu$maximum
  
  df_stu <- df_stu_opt
  
  if(result_nu$objective < log_nu[length(log_nu)]){
    df_stu <- 8
  }
  
  mod_log_n1 <- glm(formula = y ~ . , data = dat3, family = binomial(link = "logit"))
  mod_stu_n1 <- glm(formula = y ~ ., data = dat3, family = binomial(link = Gosset(df_stu)))
  step_log <- step(mod_log_n1, direction = "backward", trace = F, k=log(nrow(dat3)))
  step_stu <- step(mod_stu_n1, direction = "backward", trace = F, k=log(nrow(dat3)))
  
  step_log_b <- step(mod_log_n1, direction = "both", trace = F, k=log(nrow(dat3)))
  step_stu_b <- step(mod_stu_n1, direction = "both", trace = F, k=log(nrow(dat3)))
  
  
  mod_log_n0 <- glm(formula = y ~ 1, data = dat3, family = binomial(link = "logit"))
  mod_stu_n0 <- glm(formula = y ~ 1, data = dat3, family = binomial(link = Gosset(df_stu)))
  
  
  step_log_f <- step(mod_log_n0, direction = "forward", trace = F, k=log(nrow(dat3)), scope=(formula(mod_log_n1)))
  step_stu_f <- step(mod_stu_n0, direction = "forward", trace = F, k=log(nrow(dat3)), scope=(formula(mod_stu_n1)))
  
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
  
  step_log_c_f <- step_log_f$coefficients
  step_stu_c_f <- step_stu_f$coefficients
  pv_log <- summary(mod_log_n1)$coefficients[,4]
  
  dat_sep_f <- bind_rows(
    as.data.frame(cbind(t(pv_log),link = "BB")),
    as.data.frame(cbind(t(step_log_c_f),link = "logit")),
    as.data.frame(cbind(t(step_stu_c_f),link = "stu")))
  dat_sep_f <- dat_sep_f[-1,]
  dat_sep_f[is.na(dat_sep_f)] <- 0
  dat_b_f <- dat_sep_f
  
  step_log_c_b <- step_log_b$coefficients
  step_stu_c_b <- step_stu_b$coefficients
  pv_log <- summary(mod_log_n1)$coefficients[,4]
  
  dat_sep_b <- bind_rows(
    as.data.frame(cbind(t(pv_log),link = "BB")),
    as.data.frame(cbind(t(step_log_c_b),link = "logit")),
    as.data.frame(cbind(t(step_stu_c_b),link = "stu")))
  dat_sep_b <- dat_sep_b[-1,]
  dat_sep_b[is.na(dat_sep_b)] <- 0
  dat_b_b <- dat_sep_b
  
  print(ggarrange(plot1_po_pn,plot1_po_n,plot1_o_pn,plot1_o_n,ncol = 2,nrow = 2))
  
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
    coef_sel = dat_b,
    coef_sel_f = dat_b_f,
    coef_sel_b = dat_b_b,
    max_grille = max_grille,df_stu_opt = df_stu_opt,
    df_stu = df_stu
  ))
}



# 0.05 --------------------------------------------------------------------

nsim <- 100
C05df_stu <- C05coefs_b <- C05coefs_f <- C05coefs <- C05AICsteplog <- C05BICstepstu <- C05AICstepstu <-  C05BICsteplog <- C05max_grille <- C05df_stu_opt <- NULL
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    
    mtry <- try(num_exp_8(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable))
    variable_in <- variable + 8000
    while (inherits(mtry, "try-error")) {
      variable_in <- variable_in+1
      mtry <- try(num_exp_8(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable_in))
    } 
    
    mode1 <- mtry
    # mode1 <- num_exp_8(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable)
    C05coefs <- rbind(C05coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    C05coefs_f <- rbind(C05coefs_f,cbind(mode1[["coef_sel_f"]],rbind(outlier,outlier)))
    C05coefs_b <- rbind(C05coefs_b,cbind(mode1[["coef_sel_b"]],rbind(outlier,outlier)))
    C05AICsteplog <- rbind(C05AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    C05AICstepstu <- rbind(C05AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    C05BICsteplog <- rbind(C05BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    C05BICstepstu <- rbind(C05BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    C05max_grille <- rbind(C05max_grille,cbind(mode1[["max_grille"]],outlier))
    C05df_stu_opt <- rbind(C05df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    C05df_stu <- rbind(C05df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

cbind(C05max_grille, C05df_stu_opt, C05df_stu)
cbind(C05max_grille, C05df_stu_opt)

head(C05coefs)
C05dat5 <- C05coefs
C05dat5[C05dat5==0] <- NA
C05coefs <- C05dat5
colnames(C05coefs)[ncol(C05coefs)] <- "out"
C05summ_coef <- C05coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(C05summ_coef)

C05_SUM <- cbind(C05summ_coef[,1:2],
  exva = rowMeans(C05summ_coef[,4:7]),
  novar = rowMeans(C05summ_coef[,8:11]))

C05BIC_PLOT <- ggplot(data = C05_SUM, aes(x=out,color=link)) +
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

C05BIC_PLOT




# 0.1 --------------------------------------------------------------------

nsim <- 100
C1df_stu <- C1coefs_b <- C1coefs_f <- C1coefs <- C1AICsteplog <- C1BICstepstu <- C1AICstepstu <-  C1BICsteplog <- C1max_grille <- C1df_stu_opt <- NULL
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mtry <- try(num_exp_8(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable))
    variable_in <- variable + 8000
    while (inherits(mtry, "try-error")) {
      variable_in <- variable_in+1
      mtry <- try(num_exp_8(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable_in))
    } 
    
    mode1 <- mtry
    # mode1 <- num_exp_8(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable+20)
    C1coefs <- rbind(C1coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    C1coefs_f <- rbind(C1coefs_f,cbind(mode1[["coef_sel_f"]],rbind(outlier,outlier)))
    C1coefs_b <- rbind(C1coefs_b,cbind(mode1[["coef_sel_b"]],rbind(outlier,outlier)))
    C1AICsteplog <- rbind(C1AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    C1AICstepstu <- rbind(C1AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    C1BICsteplog <- rbind(C1BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    C1BICstepstu <- rbind(C1BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    C1max_grille <- rbind(C1max_grille,cbind(mode1[["max_grille"]],outlier))
    C1df_stu_opt <- rbind(C1df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    C1df_stu <- rbind(C1df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

cbind(C1max_grille, C1df_stu_opt)

degree_info <- data.frame(cbind(C1max_grille, C1df_stu_opt, C1df_stu[-nrow(C1df_stu)]))
colnames(degree_info) <- c("maxgri","outlier","maxalgo","out2","selected")
degree_info$sel_algo <- c(degree_info$maxalgo == degree_info$selected)
head(degree_info)

ggplot(data = degree_info, mapping = aes(x=maxgri,y=maxalgo, shape = sel_algo, colour = outlier))+
  # geom_point()+
  coord_cartesian(xlim=c(0.25,1))+
  geom_jitter()

ggplot(data = degree_info, mapping = aes(x=factor(outlier),y=selected))+
  geom_jitter(height = 0.2)
  # coord_cartesian(xlim=c(0.25,1))+
  # geom_boxplot()


C1dat5 <- C1coefs
C1dat5[C1dat5==0] <- NA
C1coefs <- C1dat5
colnames(C1coefs)[ncol(C1coefs)] <- "out"
C1summ_coef <- C1coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(C1summ_coef)

C1_SUM <- cbind(C1summ_coef[,1:2],
  exva = rowMeans(C1summ_coef[,4:7]),
  novar = rowMeans(C1summ_coef[,8:11]))

C1BIC_PLOT <- ggplot(data = C1_SUM, aes(x=out,color=link)) +
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

C1BIC_PLOT




# 0.3 --------------------------------------------------------------------

nsim <- 100
C3df_stu <- C3coefs_b <- C3coefs_f <- C3coefs <- C3AICsteplog <- C3BICstepstu <- C3AICstepstu <-  C3BICsteplog <- C3max_grille <- C3df_stu_opt <- NULL
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mtry <- try(num_exp_8(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable))
    variable_in <- variable + 8000
    while (inherits(mtry, "try-error")) {
      variable_in <- variable_in+1
      mtry <- try(num_exp_8(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable_in))
    } 
    
    mode1 <- mtry
    # mode1 <- num_exp_8(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable)
    C3coefs <- rbind(C3coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    C3coefs_f <- rbind(C3coefs_f,cbind(mode1[["coef_sel_f"]],rbind(outlier,outlier)))
    C3coefs_b <- rbind(C3coefs_b,cbind(mode1[["coef_sel_b"]],rbind(outlier,outlier)))
    C3AICsteplog <- rbind(C3AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    C3AICstepstu <- rbind(C3AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    C3BICsteplog <- rbind(C3BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    C3BICstepstu <- rbind(C3BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    C3max_grille <- rbind(C3max_grille,cbind(mode1[["max_grille"]],outlier))
    C3df_stu_opt <- rbind(C3df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    C3df_stu <- rbind(C3df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

cbind(C3max_grille, C3df_stu_opt)

C3dat5 <- C3coefs
C3dat5[C3dat5==0] <- NA
C3coefs <- C3dat5
colnames(C3coefs)[ncol(C3coefs)] <- "out"
C3summ_coef <- C3coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(C3summ_coef)

C3_SUM <- cbind(C3summ_coef[,1:2],
  exva = rowMeans(C3summ_coef[,4:7]),
  novar = rowMeans(C3summ_coef[,8:11]))

C3BIC_PLOT <- ggplot(data = C3_SUM, aes(x=out,color=link)) +
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

C3BIC_PLOT



# 0.8 --------------------------------------------------------------------

nsim <- 100
C8df_stu <- C8coefs_b <- C8coefs_f <- C8coefs <- C8AICsteplog <- C8BICstepstu <- C8AICstepstu <-  C8BICsteplog <- C8max_grille <- C8df_stu_opt <- NULL
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mtry <- try(num_exp_8(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable))
    variable_in <- variable + 8000
    while (inherits(mtry, "try-error")) {
      variable_in <- variable_in+1
      mtry <- try(num_exp_8(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable_in))
    } 
    
    mode1 <- mtry
    # mode1 <- num_exp_8(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable)
    C8coefs <- rbind(C8coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    C8coefs_f <- rbind(C8coefs_f,cbind(mode1[["coef_sel_f"]],rbind(outlier,outlier)))
    C8coefs_b <- rbind(C8coefs_b,cbind(mode1[["coef_sel_b"]],rbind(outlier,outlier)))
    C8AICsteplog <- rbind(C8AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    C8AICstepstu <- rbind(C8AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    C8BICsteplog <- rbind(C8BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    C8BICstepstu <- rbind(C8BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    C8max_grille <- rbind(C8max_grille,cbind(mode1[["max_grille"]],outlier))
    C8df_stu_opt <- rbind(C8df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    C8df_stu <- rbind(C8df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

cbind(C8max_grille, C8df_stu_opt)

C8dat5 <- C8coefs
C8dat5[C8dat5==0] <- NA
C8coefs <- C8dat5
colnames(C8coefs)[ncol(C8coefs)] <- "out"
C8summ_coef <- C8coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(C8summ_coef)

C8_SUM <- cbind(C8summ_coef[,1:2],
  exva = rowMeans(C8summ_coef[,4:7]),
  novar = rowMeans(C8summ_coef[,8:11]))

C8BIC_PLOT <- ggplot(data = C8_SUM, aes(x=out,color=link)) +
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

C8BIC_PLOT


ggarrange(C05BIC_PLOT, C1BIC_PLOT, C3BIC_PLOT, C8BIC_PLOT, ncol = 2, nrow = 2) 
