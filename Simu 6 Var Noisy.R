library(ggplot2)
# Numerical illustration --------------------------------------------------

plot1_po_pn <- plot1_po_n <- plot1_o_pn <- plot1_o_n <- NULL
set.seed(30) # 4  ou 11
g_deviance <- function(nu, data_in) {
  logLik(glm(y~.,family = binomial(Gosset(nu)),data = data_in, maxit = 100))
}

num_exp6 <- function(n_obs,s,n_out,df_stu,seeed,out_y="meth2"){
  set.seed(seeed)
  
  n_obs <- n_obs*1.3
  
  alpha <- 0
  beta1 <- c(.8, 0.4, -0.4, 0.2, rep(0, 2))
  Sigma <- matrix(0, nrow = 6, ncol = 6)
  diag(Sigma) <- 1
  Ex <- cbind(MASS::mvrnorm(n_obs, rep(0, 6), Sigma))
  eta1 <- alpha + as.matrix(Ex) %*% (beta1)
  h_x_beta <- S(eta1,s)
  y <- rbinom(n = n_obs,size = 1,prob = h_x_beta)  
  dat2 <- data.frame(y,Ex)
  dat2$y <- as.factor(dat2$y)
  dat3 <- dat2
  
  n_out_f <- n_out
  Sigma2 <- matrix(0, nrow = 6, ncol = 6)
  diag(Sigma2) <- 0.5
  if(n_out_f > 0){
    if(n_out_f==1){
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = 2, c(1.6,0.8,-0.8,0.4,2,2), Sigma2)),y=rep(0,n_out_f)) 
      dat3 <- rbind(dat2,outliers[1,])
    }else{
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = n_out_f, c(1.6,0.8,-0.8,0.4,2,2), Sigma2)),y=rep(0,n_out_f)) 
      dat3 <- rbind(dat2,outliers)
    }
  }
  
  dat_test <- dat3[1:30,]
  n_test <- nrow(dat_test)
  dat3 <- dat3[31:nrow(dat3),]
  
  nu_expl <- c(1,8)
  log_nu <- NULL
  i <- 1
  log_nu1 <- tryCatch(
    for (nu in nu_expl) {
      log_nu[i] <- logLik(glm(formula = y ~. , data = dat3, family = binomial(link = Gosset(nu)),maxit = 1000))
      i <- i+1
    },
    warning = 
      function(e){
        return(0)
      })
  # No models without convergence
  if(length(log_nu1) == 1 ){stop("")}
  
  log_nu <- NULL
  i <- 1
  for (nu in nu_expl) {
    log_nu[i] <- logLik(glm(formula = y ~. , data = dat3, family = binomial(link = Gosset(nu)),maxit = 1000))
    i <- i+1
  }
  
  log_1 <- log_nu[nu_expl==1]
  log_8 <- log_nu[nu_expl==8]
  
  df_stu <- 8
  
  if(log_1 > log_8){
    result_nu <- optimize(g_deviance, c(0.25, 1), maximum = T, dat3)
    df_stu_opt <- result_nu$maximum
    df_stu <- df_stu_opt
  }
  
  mod_log_n1 <- glm(formula = y ~ . , data = dat3, family = binomial(link = "logit"))
  
  log_nu3 <- tryCatch(
    glm(formula = y ~ ., data = dat3, family = binomial(link = Gosset(df_stu))),
    warning = 
      function(e){
        return(0)
      })
  # No models without convergence
  if(length(log_nu3) == 1 ){stop("")}
  
  mod_stu_n1 <- glm(formula = y ~ ., data = dat3, family = binomial(link = Gosset(df_stu)))
  step_log <- step(mod_log_n1, direction = "backward", trace = F, k=log(nrow(dat3)))
  step_stu <- step(mod_stu_n1, direction = "backward", trace = F, k=log(nrow(dat3)))
  
  pred_log <- ifelse(predict(step_log, dat_test, type = "response")<0.5,0,1)
  tab_log <- table(dat_test$y,pred_log)
  acc_log <- sum(diag(tab_log))/n_test
  
  pred_stu <- ifelse(predict(step_stu, dat_test, type = "response")<0.5,0,1)
  tab_stu <- table(dat_test$y,pred_stu)
  acc_stu <- sum(diag(tab_stu))/n_test
  
  
  # step_log_b <- step(mod_log_n1, direction = "both", trace = F, k=log(nrow(dat3)))
  # step_stu_b <- step(mod_stu_n1, direction = "both", trace = F, k=log(nrow(dat3)))
  
  # mod_log_n0 <- glm(formula = y ~ 1, data = dat3, family = binomial(link = "logit"))
  # mod_stu_n0 <- glm(formula = y ~ 1, data = dat3, family = binomial(link = Gosset(df_stu)))
  
  # step_log_f <- step(mod_log_n0, direction = "forward", trace = F, k=log(nrow(dat3)), scope=(formula(mod_log_n1)))
  # step_stu_f <- step(mod_stu_n0, direction = "forward", trace = F, k=log(nrow(dat3)), scope=(formula(mod_stu_n1)))
  
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
  
  # step_log_c_b <- step_log_b$coefficients
  # step_stu_c_b <- step_stu_b$coefficients
  # pv_log <- summary(mod_log_n1)$coefficients[,4]
  # 
  # dat_sep_b <- bind_rows(
  #   as.data.frame(cbind(t(pv_log),link = "BB")),
  #   as.data.frame(cbind(t(step_log_c_b),link = "logit")),
  #   as.data.frame(cbind(t(step_stu_c_b),link = "stu")))
  # dat_sep_b <- dat_sep_b[-1,]
  # dat_sep_b[is.na(dat_sep_b)] <- 0
  # dat_b_b <- dat_sep_b
  
  
  
  return(list(
    acc_log = acc_log,
    acc_stu = acc_stu,
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
    # coef_sel_b = dat_b_b,
    # max_grille = max_grille,
    df_stu_opt = df_stu_opt,
    df_stu = df_stu
  ))
}



# 0.05 --------------------------------------------------------------------

nsim <- 100
C005acc_stu <- C005acc_log <- C005df_stu <- C005coefs_b <- C005coefs_f <- C005coefs <- C005AICsteplog <- C005BICstepstu <- C005AICstepstu <-  C005BICsteplog <- C005max_grille <- C005df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp6(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp6(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp6(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable)
    C005coefs <- rbind(C005coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    C005acc_log <- rbind(C005acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    C005acc_stu <- rbind(C005acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    C005AICsteplog <- rbind(C005AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    C005AICstepstu <- rbind(C005AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    C005BICsteplog <- rbind(C005BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    C005BICstepstu <- rbind(C005BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # C005df_stu_opt <- rbind(C005df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    C005df_stu <- rbind(C005df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(C005coefs)
dat5 <- C005coefs[,-ncol(C005coefs)]
dat5[dat5==0] <- NA
C005coefs[,-ncol(C005coefs)] <- dat5
colnames(C005coefs)[ncol(C005coefs)] <- "out"

C005summ_coef <- C005coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(C005summ_coef)

save(C005summ_coef, file = "data/Simu 6 Var Noisy/C005summ_coef.RData")

C005_SUM <- cbind(C005summ_coef[,1:2],
  exva = rowMeans(C005summ_coef[,4:7]),
  novar = rowMeans(C005summ_coef[,8:9]))

C005summ_coef_p <- ggplot(data = C005_SUM, aes(x=out,color=link)) +
  geom_line(aes(y=novar, linetype = "dashed")) +
  geom_line(aes(y=exva,linetype = "solid")) +
  geom_point(aes(y=novar))+
  geom_point(aes(y=exva))+
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.05$")

C005summ_coef_p

C005acc_logit <- as.data.frame(C005acc_log)
colnames(C005acc_logit)[1] <- "loglik"
C005acc_logit$link <- "logit"

C005acc_stu <- as.data.frame(C005acc_stu)
colnames(C005acc_stu)[1] <- "loglik"
C005acc_stu$link <- "stu"

C005accliks <- rbind(C005acc_logit,C005acc_stu)
head(C005accliks)
colnames(C005accliks) <- c("loglik","outlier","link")

save(C005accliks, file = "data/Simu 6 Var Noisy/C005accliks.RData")

C005acc_plot <- ggplot(data = C005accliks[C005accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.05")

C005acc_plot

# 0.1 --------------------------------------------------------------------

nsim <- 100
C01acc_stu <- C01acc_log <- C01df_stu <- C01coefs_b <- C01coefs_f <- C01coefs <- C01AICsteplog <- C01BICstepstu <- C01AICstepstu <-  C01BICsteplog <- C01max_grille <- C01df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp6(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp6(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp6(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable)
    C01coefs <- rbind(C01coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    C01acc_log <- rbind(C01acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    C01acc_stu <- rbind(C01acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    C01AICsteplog <- rbind(C01AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    C01AICstepstu <- rbind(C01AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    C01BICsteplog <- rbind(C01BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    C01BICstepstu <- rbind(C01BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # C01df_stu_opt <- rbind(C01df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    C01df_stu <- rbind(C01df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(C01coefs)
dat5 <- C01coefs[,-ncol(C01coefs)]
dat5[dat5==0] <- NA
C01coefs[,-ncol(C01coefs)] <- dat5
colnames(C01coefs)[ncol(C01coefs)] <- "out"

C01summ_coef <- C01coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(C01summ_coef)

save(C01summ_coef, file = "data/Simu 6 Var Noisy/C01summ_coef.RData")

C01_SUM <- cbind(C01summ_coef[,1:2],
  exva = rowMeans(C01summ_coef[,4:7]),
  novar = rowMeans(C01summ_coef[,8:9]))

C01summ_coef_p <- ggplot(data = C01_SUM, aes(x=out,color=link)) +
  geom_line(aes(y=novar, linetype = "dashed")) +
  geom_line(aes(y=exva,linetype = "solid")) +
  geom_point(aes(y=novar))+
  geom_point(aes(y=exva))+
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.1$")

C01summ_coef_p

C01acc_logit <- as.data.frame(C01acc_log)
colnames(C01acc_logit)[1] <- "loglik"
C01acc_logit$link <- "logit"

C01acc_stu <- as.data.frame(C01acc_stu)
colnames(C01acc_stu)[1] <- "loglik"
C01acc_stu$link <- "stu"

C01accliks <- rbind(C01acc_logit,C01acc_stu)
head(C01accliks)
colnames(C01accliks) <- c("loglik","outlier","link")

save(C01accliks, file = "data/Simu 6 Var Noisy/C01accliks.RData")

C01acc_plot <- ggplot(data = C01accliks[C01accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.1")

C01acc_plot

# 0.3 --------------------------------------------------------------------
nsim <- 100
C03acc_stu <- C03acc_log <- C03df_stu <- C03coefs_b <- C03coefs_f <- C03coefs <- C03AICsteplog <- C03BICstepstu <- C03AICstepstu <-  C03BICsteplog <- C03max_grille <- C03df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp6(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp6(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp6(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable)
    C03coefs <- rbind(C03coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    C03acc_log <- rbind(C03acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    C03acc_stu <- rbind(C03acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    C03AICsteplog <- rbind(C03AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    C03AICstepstu <- rbind(C03AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    C03BICsteplog <- rbind(C03BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    C03BICstepstu <- rbind(C03BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # C03df_stu_opt <- rbind(C03df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    C03df_stu <- rbind(C03df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(C03coefs)
dat5 <- C03coefs[,-ncol(C03coefs)]
dat5[dat5==0] <- NA
C03coefs[,-ncol(C03coefs)] <- dat5
colnames(C03coefs)[ncol(C03coefs)] <- "out"

C03summ_coef <- C03coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(C03summ_coef)

save(C03summ_coef, file = "data/Simu 6 Var Noisy/C03summ_coef.RData")

C03_SUM <- cbind(C03summ_coef[,1:2],
  exva = rowMeans(C03summ_coef[,4:7]),
  novar = rowMeans(C03summ_coef[,8:9]))

C03summ_coef_p <- ggplot(data = C03_SUM, aes(x=out,color=link)) +
  geom_line(aes(y=novar, linetype = "dashed")) +
  geom_line(aes(y=exva,linetype = "solid")) +
  geom_point(aes(y=novar))+
  geom_point(aes(y=exva))+
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.3$")

C03summ_coef_p

C03acc_logit <- as.data.frame(C03acc_log)
colnames(C03acc_logit)[1] <- "loglik"
C03acc_logit$link <- "logit"

C03acc_stu <- as.data.frame(C03acc_stu)
colnames(C03acc_stu)[1] <- "loglik"
C03acc_stu$link <- "stu"

C03accliks <- rbind(C03acc_logit,C03acc_stu)
head(C03accliks)
colnames(C03accliks) <- c("loglik","outlier","link")

save(C03accliks, file = "data/Simu 6 Var Noisy/C03accliks.RData")

C03acc_plot <- ggplot(data = C03accliks[C03accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.3")

C03acc_plot

# 0.8 --------------------------------------------------------------------

nsim <- 100
C08acc_stu <- C08acc_log <- C08df_stu <- C08coefs_b <- C08coefs_f <- C08coefs <- C08AICsteplog <- C08BICstepstu <- C08AICstepstu <-  C08BICsteplog <- C08max_grille <- C08df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp6(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp6(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp6(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable)
    C08coefs <- rbind(C08coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    C08acc_log <- rbind(C08acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    C08acc_stu <- rbind(C08acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    C08AICsteplog <- rbind(C08AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    C08AICstepstu <- rbind(C08AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    C08BICsteplog <- rbind(C08BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    C08BICstepstu <- rbind(C08BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # C08df_stu_opt <- rbind(C08df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    C08df_stu <- rbind(C08df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(C08coefs)
dat5 <- C08coefs[,-ncol(C08coefs)]
dat5[dat5==0] <- NA
C08coefs[,-ncol(C08coefs)] <- dat5
colnames(C08coefs)[ncol(C08coefs)] <- "out"

C08summ_coef <- C08coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(C08summ_coef)

save(C08summ_coef, file = "data/Simu 6 Var Noisy/C08summ_coef.RData")

C08_SUM <- cbind(C08summ_coef[,1:2],
  exva = rowMeans(C08summ_coef[,4:7]),
  novar = rowMeans(C08summ_coef[,8:9]))

C08summ_coef_p <- ggplot(data = C08_SUM, aes(x=out,color=link)) +
  geom_line(aes(y=novar, linetype = "dashed")) +
  geom_line(aes(y=exva,linetype = "solid")) +
  geom_point(aes(y=novar))+
  geom_point(aes(y=exva))+
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.8$")

C08summ_coef_p

C08acc_logit <- as.data.frame(C08acc_log)
colnames(C08acc_logit)[1] <- "loglik"
C08acc_logit$link <- "logit"

C08acc_stu <- as.data.frame(C08acc_stu)
colnames(C08acc_stu)[1] <- "loglik"
C08acc_stu$link <- "stu"

C08accliks <- rbind(C08acc_logit,C08acc_stu)
head(C08accliks)
colnames(C08accliks) <- c("loglik","outlier","link")

save(C08accliks, file = "data/Simu 6 Var Noisy/C08accliks.RData")

C08acc_plot <- ggplot(data = C08accliks[C08accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.8")

C08acc_plot


# Plots -------------------------------------------------------------------


tikz(file = "tex/Av_coef_6_Noisy.tex", standAlone=F, width = 7.5, height = 4.5)
# endoffile <- dev.off()

ggarrange(C005summ_coef_p, C01summ_coef_p, C03summ_coef_p, C08summ_coef_p, ncol = 2, nrow = 2) 
ggarrange(C01acc_plot, C01acc_plot, C03acc_plot, C08acc_plot, ncol = 2, nrow = 2) 
