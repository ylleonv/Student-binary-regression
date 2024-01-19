library(ggplot2)
# Numerical illustration --------------------------------------------------

plot1_po_pn <- plot1_po_n <- plot1_o_pn <- plot1_o_n <- NULL
set.seed(30) # 4  ou 11
g_deviance <- function(nu, data_in) {
  logLik(glm(y~.,family = binomial(Gosset(nu)),data = data_in, maxit = 100))
}

num_exp_5 <- function(n_obs,s,n_out,df_stu,seeed,out_y="meth2"){
  set.seed(seeed)
  alpha <- 0
  # beta1 <- c(1, rep(0, 1))
  Sigma <- matrix(0, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  beta1 <- c(.8, -.6, rep(0, 2))
  
  Sigma <- matrix(0, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  Ex <- cbind(MASS::mvrnorm(n_obs, rep(0, 4), Sigma))
  
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
  Sigma2 <- matrix(0, nrow = 4, ncol = 4)
  diag(Sigma2) <- 0.25
  if(n_out_f > 0){
    if(n_out_f==1){
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = 2, c(1,-1,2,2), Sigma2)),y=rep(0,n_out_f)) 
      dat3 <- rbind(dat2,outliers[1,])
    }else{
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = n_out_f, c(1,-1,2,2), Sigma2)),y=rep(0,n_out_f)) 
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
  
  # print(result_nu$objective)
  # print(log_nu[length(log_nu)])
  # 
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
  
  
  # step_log_f <- step(mod_log_n0, direction = "forward", trace = F, k=log(nrow(dat3)), scope=(formula(mod_log_n1)))
  # step_stu_f <- step(mod_stu_n0, direction = "forward", trace = F, k=log(nrow(dat3)), scope=(formula(mod_stu_n1)))
  # 
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
  
  # step_log_c_f <- step_log_f$coefficients
  # step_stu_c_f <- step_stu_f$coefficients
  # pv_log <- summary(mod_log_n1)$coefficients[,4]
  # 
  # dat_sep_f <- bind_rows(
  #   as.data.frame(cbind(t(pv_log),link = "BB")),
  #   as.data.frame(cbind(t(step_log_c_f),link = "logit")),
  #   as.data.frame(cbind(t(step_stu_c_f),link = "stu")))
  # dat_sep_f <- dat_sep_f[-1,]
  # dat_sep_f[is.na(dat_sep_f)] <- 0
  # dat_b_f <- dat_sep_f
  
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
    coef_sel = dat_b,
    # coef_sel_f = dat_b_f,
    coef_sel_b = dat_b_b,
    max_grille = max_grille,
    df_stu_opt = df_stu_opt,
    df_stu = df_stu
  ))
}



# 0.05 --------------------------------------------------------------------

nsim <- 100
B05df_stu <- B05coefs_b <- B05coefs_f <- B05coefs <- B05AICsteplog <- B05BICstepstu <- B05AICstepstu <-  B05BICsteplog <- B05max_grille <- B05df_stu_opt <- NULL
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    
    mtry <- try(num_exp_5(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable))
    variable_in <- variable + 8000
    while (inherits(mtry, "try-error")) {
      variable_in <- variable_in+1
      mtry <- try(num_exp_5(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable_in))
    } 
    
    mode1 <- mtry
    # mode1 <- num_exp_5(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable)
    B05coefs <- rbind(B05coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    # B05coefs_f <- rbind(B05coefs_f,cbind(mode1[["coef_sel_f"]],rbind(outlier,outlier)))
    B05coefs_b <- rbind(B05coefs_b,cbind(mode1[["coef_sel_b"]],rbind(outlier,outlier)))
    B05AICsteplog <- rbind(B05AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    B05AICstepstu <- rbind(B05AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    B05BICsteplog <- rbind(B05BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    B05BICstepstu <- rbind(B05BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    B05max_grille <- rbind(B05max_grille,cbind(mode1[["max_grille"]],outlier))
    B05df_stu_opt <- rbind(B05df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    B05df_stu <- rbind(B05df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

view(B05df_stu)
cbind(B05max_grille, B05df_stu_opt, B05df_stu)

head(B05coefs)
B05dat5 <- B05coefs
B05dat5[B05dat5==0] <- NA
B05coefs <- B05dat5
colnames(B05coefs)[ncol(B05coefs)] <- "out"
B05summ_coef <- B05coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(B05summ_coef)

B05_SUM <- cbind(B05summ_coef[,1:2],
  exva = rowMeans(B05summ_coef[,4:5]),
  novar = rowMeans(B05summ_coef[,6:7]))

B05BIC_PLOT <- ggplot(data = B05_SUM, aes(x=out,color=link)) +
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

B05BIC_PLOT




# 0.1 --------------------------------------------------------------------

nsim <- 100
B1df_stu <- B1coefs_b <- B1coefs_f <- B1coefs <- B1AICsteplog <- B1BICstepstu <- B1AICstepstu <-  B1BICsteplog <- B1max_grille <- B1df_stu_opt <- NULL
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mtry <- try(num_exp_5(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable))
    variable_in <- variable + 8000
    while (inherits(mtry, "try-error")) {
      variable_in <- variable_in+1
      mtry <- try(num_exp_5(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable_in))
    } 
    
    mode1 <- mtry
    # mode1 <- num_exp_5(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable+20)
    B1coefs <- rbind(B1coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    # B1coefs_f <- rbind(B1coefs_f,cbind(mode1[["coef_sel_f"]],rbind(outlier,outlier)))
    B1coefs_b <- rbind(B1coefs_b,cbind(mode1[["coef_sel_b"]],rbind(outlier,outlier)))
    B1AICsteplog <- rbind(B1AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    B1AICstepstu <- rbind(B1AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    B1BICsteplog <- rbind(B1BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    B1BICstepstu <- rbind(B1BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    B1max_grille <- rbind(B1max_grille,cbind(mode1[["max_grille"]],outlier))
    B1df_stu_opt <- rbind(B1df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    B1df_stu <- rbind(B1df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}



View(B1df_stu_opt)
degree_info <- data.frame(cbind(B1max_grille, B1df_stu_opt, B1df_stu[-nrow(B1df_stu)]))
colnames(degree_info) <- c("maxgri","outlier","maxalgo","out2","selected")




B1dat5 <- B1coefs
B1dat5[B1dat5==0] <- NA
B1coefs <- B1dat5
colnames(B1coefs)[ncol(B1coefs)] <- "out"
B1summ_coef <- B1coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(B1summ_coef)

B1_SUM <- cbind(B1summ_coef[,1:2],
  exva = rowMeans(B1summ_coef[,4:5]),
  novar = rowMeans(B1summ_coef[,6:7]))

B1BIC_PLOT <- ggplot(data = B1_SUM, aes(x=out,color=link)) +
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

B1BIC_PLOT




# 0.3 --------------------------------------------------------------------

nsim <- 100
B3df_stu <- B3coefs_b <- B3coefs_f <- B3coefs <- B3AICsteplog <- B3BICstepstu <- B3AICstepstu <-  B3BICsteplog <- B3max_grille <- B3df_stu_opt <- NULL
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mtry <- try(num_exp_5(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable))
    variable_in <- variable + 8000
    while (inherits(mtry, "try-error")) {
      variable_in <- variable_in+1
      mtry <- try(num_exp_5(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable_in))
    } 
    
    mode1 <- mtry
    # mode1 <- num_exp_5(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable)
    B3coefs <- rbind(B3coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    # B3coefs_f <- rbind(B3coefs_f,cbind(mode1[["coef_sel_f"]],rbind(outlier,outlier)))
    B3coefs_b <- rbind(B3coefs_b,cbind(mode1[["coef_sel_b"]],rbind(outlier,outlier)))
    B3AICsteplog <- rbind(B3AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    B3AICstepstu <- rbind(B3AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    B3BICsteplog <- rbind(B3BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    B3BICstepstu <- rbind(B3BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    B3max_grille <- rbind(B3max_grille,cbind(mode1[["max_grille"]],outlier))
    B3df_stu_opt <- rbind(B3df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    B3df_stu <- rbind(B3df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

cbind(B3max_grille, B3df_stu_opt)

degree_info <- data.frame(cbind(B3max_grille, B3df_stu_opt, B3df_stu[-nrow(B3df_stu)]))
colnames(degree_info) <- c("maxgri","outlier","maxalgo","out2","selected")

B3dat5 <- B3coefs
B3dat5[B3dat5==0] <- NA
B3coefs <- B3dat5
colnames(B3coefs)[ncol(B3coefs)] <- "out"
B3summ_coef <- B3coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(B3summ_coef)

B3_SUM <- cbind(B3summ_coef[,1:2],
  exva = rowMeans(B3summ_coef[,4:5]),
  novar = rowMeans(B3summ_coef[,6:7]))

B3BIC_PLOT <- ggplot(data = B3_SUM, aes(x=out,color=link)) +
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

B3BIC_PLOT



# 0.8 --------------------------------------------------------------------

nsim <- 100
B8df_stu <- B8coefs_b <- B8coefs_f <- B8coefs <- B8AICsteplog <- B8BICstepstu <- B8AICstepstu <-  B8BICsteplog <- B8max_grille <- B8df_stu_opt <- NULL
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mtry <- try(num_exp_5(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable))
    variable_in <- variable + 8000
    while (inherits(mtry, "try-error")) {
      variable_in <- variable_in+1
      mtry <- try(num_exp_5(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable_in))
    } 
    
    mode1 <- mtry
    # mode1 <- num_exp_5(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable)
    B8coefs <- rbind(B8coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    # B8coefs_f <- rbind(B8coefs_f,cbind(mode1[["coef_sel_f"]],rbind(outlier,outlier)))
    B8coefs_b <- rbind(B8coefs_b,cbind(mode1[["coef_sel_b"]],rbind(outlier,outlier)))
    B8AICsteplog <- rbind(B8AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    B8AICstepstu <- rbind(B8AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    B8BICsteplog <- rbind(B8BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    B8BICstepstu <- rbind(B8BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    B8max_grille <- rbind(B8max_grille,cbind(mode1[["max_grille"]],outlier))
    B8df_stu_opt <- rbind(B8df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    B8df_stu <- rbind(B8df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

cbind(B8max_grille, B8df_stu_opt)

degree_info <- data.frame(cbind(B8max_grille, B8df_stu_opt, B8df_stu[-nrow(B8df_stu)]))
colnames(degree_info) <- c("maxgri","outlier","maxalgo","out2","selected")

B8dat5 <- B8coefs
B8dat5[B8dat5==0] <- NA
B8coefs <- B8dat5
colnames(B8coefs)[ncol(B8coefs)] <- "out"
B8summ_coef <- B8coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(B8summ_coef)

B8_SUM <- cbind(B8summ_coef[,1:2],
  exva = rowMeans(B8summ_coef[,4:5]),
  novar = rowMeans(B8summ_coef[,6:7]))

B8BIC_PLOT <- ggplot(data = B8_SUM, aes(x=out,color=link)) +
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

B8BIC_PLOT


# Graficos ----------------------------------------------------------------

B05dat5 <- B05coefs_f
B05dat5[B05dat5==0] <- NA
B05coefs_f <- B05dat5
colnames(B05coefs_f)[ncol(B05coefs_f)] <- "out"
B05summ_coef <- B05coefs_f %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(B05summ_coef)

B05_SUM <- cbind(B05summ_coef[,1:2],
  exva = rowMeans(B05summ_coef[,4:5]),
  novar = rowMeans(B05summ_coef[,6:7]))

B05BIC_PLOT <- ggplot(data = B05_SUM, aes(x=out,color=link)) +
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

B05BIC_PLOT

B1dat5 <- B1coefs_f
B1dat5[B1dat5==0] <- NA
B1coefs_f <- B1dat5
colnames(B1coefs_f)[ncol(B1coefs_f)] <- "out"
B1summ_coef <- B1coefs_f %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(B1summ_coef)

B1_SUM <- cbind(B1summ_coef[,1:2],
  exva = rowMeans(B1summ_coef[,4:5]),
  novar = rowMeans(B1summ_coef[,6:7]))

B1BIC_PLOT <- ggplot(data = B1_SUM, aes(x=out,color=link)) +
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

B1BIC_PLOT

B3dat5 <- B3coefs_f
B3dat5[B3dat5==0] <- NA
B3coefs_f <- B3dat5
colnames(B3coefs_f)[ncol(B3coefs_f)] <- "out"
B3summ_coef <- B3coefs_f %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(B3summ_coef)

B3_SUM <- cbind(B3summ_coef[,1:2],
  exva = rowMeans(B3summ_coef[,4:5]),
  novar = rowMeans(B3summ_coef[,6:7]))

B3BIC_PLOT <- ggplot(data = B3_SUM, aes(x=out,color=link)) +
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

B3BIC_PLOT

B8dat5 <- B8coefs_f
B8dat5[B8dat5==0] <- NA
B8coefs_f <- B8dat5
colnames(B8coefs_f)[ncol(B8coefs_f)] <- "out"
B8summ_coef <- B8coefs_f %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(B8summ_coef)

B8_SUM <- cbind(B8summ_coef[,1:2],
  exva = rowMeans(B8summ_coef[,4:5]),
  novar = rowMeans(B8summ_coef[,6:7]))

B8BIC_PLOT <- ggplot(data = B8_SUM, aes(x=out,color=link)) +
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

B8BIC_PLOT
ggarrange(B05BIC_PLOT, B1BIC_PLOT, B3BIC_PLOT, B8BIC_PLOT, ncol = 2, nrow = 2) 
