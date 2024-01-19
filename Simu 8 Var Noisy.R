library(ggplot2)
# Numerical illustration --------------------------------------------------

plot1_po_pn <- plot1_po_n <- plot1_o_pn <- plot1_o_n <- NULL
set.seed(30) # 4  ou 11
g_deviance <- function(nu, data_in) {
  logLik(glm(y~.,family = binomial(Gosset(nu)),data = data_in, maxit = 100))
}

num_exp8 <- function(n_obs,s,n_out,df_stu,seeed,out_y="meth2"){
  set.seed(seeed)
  
  n_obs <- n_obs*1.3
  
  alpha <- 0
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
  
  n_out_f <- n_out
  Sigma2 <- matrix(0, nrow = 8, ncol = 8)
  diag(Sigma2) <- 0.5
  if(n_out_f > 0){
    if(n_out_f==1){
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = 2, c(1.6,0.8,-0.8,0.4,2,2,2,2), Sigma2)),y=rep(0,n_out_f)) 
      dat3 <- rbind(dat2,outliers[1,])
    }else{
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = n_out_f, c(1.6,0.8,-0.8,0.4,2,2,2,2), Sigma2)),y=rep(0,n_out_f)) 
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
D005acc_stu <- D005acc_log <- D005df_stu <- D005coefs_b <- D005coefs_f <- D005coefs <- D005AICsteplog <- D005BICstepstu <- D005AICstepstu <-  D005BICsteplog <- D005max_grille <- D005df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    # variable <- variable+1
    # variable_e <- variable_e+1
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp8(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp8(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable),
        silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp8(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable)
    D005coefs <- rbind(D005coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    D005acc_log <- rbind(D005acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    D005acc_stu <- rbind(D005acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    D005AICsteplog <- rbind(D005AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    D005AICstepstu <- rbind(D005AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    D005BICsteplog <- rbind(D005BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    D005BICstepstu <- rbind(D005BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # D005df_stu_opt <- rbind(D005df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    D005df_stu <- rbind(D005df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(D005coefs)
dat5 <- D005coefs[,-ncol(D005coefs)]
dat5[dat5==0] <- NA
D005coefs[,-ncol(D005coefs)] <- dat5
colnames(D005coefs)[ncol(D005coefs)] <- "out"

D005summ_coef <- D005coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(D005summ_coef)

D005_SUM <- cbind(D005summ_coef[,1:2],
  exva = rowMeans(D005summ_coef[,4:7]),
  novar = rowMeans(D005summ_coef[,8:11]))

D005summ_coef_p <- ggplot(data = D005_SUM, aes(x=out,color=link)) +
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

D005summ_coef_p

D005acc_logit <- as.data.frame(D005acc_log)
colnames(D005acc_logit)[1] <- "loglik"
D005acc_logit$link <- "logit"

D005acc_stu <- as.data.frame(D005acc_stu)
colnames(D005acc_stu)[1] <- "loglik"
D005acc_stu$link <- "stu"

D005accliks <- rbind(D005acc_logit,D005acc_stu)
head(D005accliks)
colnames(D005accliks) <- c("loglik","outlier","link")

D005acc_plot <- ggplot(data = D005accliks[D005accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.05")

D005acc_plot

# 0.1 --------------------------------------------------------------------

nsim <- 100
D01acc_stu <- D01acc_log <- D01df_stu <- D01coefs_b <- D01coefs_f <- D01coefs <- D01AICsteplog <- D01BICstepstu <- D01AICstepstu <-  D01BICsteplog <- D01max_grille <- D01df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    # variable <- variable+1
    # variable_e <- variable_e+1
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp8(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp8(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp8(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable)
    D01coefs <- rbind(D01coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    D01acc_log <- rbind(D01acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    D01acc_stu <- rbind(D01acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    D01AICsteplog <- rbind(D01AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    D01AICstepstu <- rbind(D01AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    D01BICsteplog <- rbind(D01BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    D01BICstepstu <- rbind(D01BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # D01df_stu_opt <- rbind(D01df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    D01df_stu <- rbind(D01df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(D01coefs)
dat5 <- D01coefs[,-ncol(D01coefs)]
dat5[dat5==0] <- NA
D01coefs[,-ncol(D01coefs)] <- dat5
colnames(D01coefs)[ncol(D01coefs)] <- "out"

D01summ_coef <- D01coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(D01summ_coef)

D01_SUM <- cbind(D01summ_coef[,1:2],
  exva = rowMeans(D01summ_coef[,4:7]),
  novar = rowMeans(D01summ_coef[,8:11]))

D01summ_coef_p <- ggplot(data = D01_SUM, aes(x=out,color=link)) +
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

D01summ_coef_p

D01acc_logit <- as.data.frame(D01acc_log)
colnames(D01acc_logit)[1] <- "loglik"
D01acc_logit$link <- "logit"

D01acc_stu <- as.data.frame(D01acc_stu)
colnames(D01acc_stu)[1] <- "loglik"
D01acc_stu$link <- "stu"

D01accliks <- rbind(D01acc_logit,D01acc_stu)
head(D01accliks)
colnames(D01accliks) <- c("loglik","outlier","link")

D01acc_plot <- ggplot(data = D01accliks[D01accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.1")

D01acc_plot

# 0.3 --------------------------------------------------------------------
nsim <- 100
D03acc_stu <- D03acc_log <- D03df_stu <- D03coefs_b <- D03coefs_f <- D03coefs <- D03AICsteplog <- D03BICstepstu <- D03AICstepstu <-  D03BICsteplog <- D03max_grille <- D03df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    # variable <- variable+1
    # variable_e <- variable_e+1
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp8(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp8(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp8(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable)
    D03coefs <- rbind(D03coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    D03acc_log <- rbind(D03acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    D03acc_stu <- rbind(D03acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    D03AICsteplog <- rbind(D03AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    D03AICstepstu <- rbind(D03AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    D03BICsteplog <- rbind(D03BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    D03BICstepstu <- rbind(D03BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # D03df_stu_opt <- rbind(D03df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    D03df_stu <- rbind(D03df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(D03coefs)
dat5 <- D03coefs[,-ncol(D03coefs)]
dat5[dat5==0] <- NA
D03coefs[,-ncol(D03coefs)] <- dat5
colnames(D03coefs)[ncol(D03coefs)] <- "out"

D03summ_coef <- D03coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(D03summ_coef)

D03_SUM <- cbind(D03summ_coef[,1:2],
  exva = rowMeans(D03summ_coef[,4:7]),
  novar = rowMeans(D03summ_coef[,8:11]))

D03summ_coef_p <- ggplot(data = D03_SUM, aes(x=out,color=link)) +
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

D03summ_coef_p

D03acc_logit <- as.data.frame(D03acc_log)
colnames(D03acc_logit)[1] <- "loglik"
D03acc_logit$link <- "logit"

D03acc_stu <- as.data.frame(D03acc_stu)
colnames(D03acc_stu)[1] <- "loglik"
D03acc_stu$link <- "stu"

D03accliks <- rbind(D03acc_logit,D03acc_stu)
head(D03accliks)
colnames(D03accliks) <- c("loglik","outlier","link")

D03acc_plot <- ggplot(data = D03accliks[D03accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.3")

D03acc_plot

# 0.8 --------------------------------------------------------------------

nsim <- 100
D08acc_stu <- D08acc_log <- D08df_stu <- D08coefs_b <- D08coefs_f <- D08coefs <- D08AICsteplog <- D08BICstepstu <- D08AICstepstu <-  D08BICsteplog <- D08max_grille <- D08df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    # variable <- variable+1
    # variable_e <- variable_e+1
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp8(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp8(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp8(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable)
    D08coefs <- rbind(D08coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    D08acc_log <- rbind(D08acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    D08acc_stu <- rbind(D08acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    D08AICsteplog <- rbind(D08AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    D08AICstepstu <- rbind(D08AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    D08BICsteplog <- rbind(D08BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    D08BICstepstu <- rbind(D08BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # D08df_stu_opt <- rbind(D08df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    D08df_stu <- rbind(D08df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(D08coefs)
dat5 <- D08coefs[,-ncol(D08coefs)]
dat5[dat5==0] <- NA
D08coefs[,-ncol(D08coefs)] <- dat5
colnames(D08coefs)[ncol(D08coefs)] <- "out"

D08summ_coef <- D08coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(D08summ_coef)

D08_SUM <- cbind(D08summ_coef[,1:2],
  exva = rowMeans(D08summ_coef[,4:7]),
  novar = rowMeans(D08summ_coef[,8:11]))

D08summ_coef_p <- ggplot(data = D08_SUM, aes(x=out,color=link)) +
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

D08summ_coef_p

D08acc_logit <- as.data.frame(D08acc_log)
colnames(D08acc_logit)[1] <- "loglik"
D08acc_logit$link <- "logit"

D08acc_stu <- as.data.frame(D08acc_stu)
colnames(D08acc_stu)[1] <- "loglik"
D08acc_stu$link <- "stu"

D08accliks <- rbind(D08acc_logit,D08acc_stu)
head(D08accliks)
colnames(D08accliks) <- c("loglik","outlier","link")

D08acc_plot <- ggplot(data = D08accliks[D08accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.8")

D08acc_plot


# Plots -------------------------------------------------------------------

tikz(file = "tex/Av_coef_8_Noisy.tex", standAlone=F, width = 7.5, height = 4.5)
# endoffile <- dev.off()

ggarrange(D005summ_coef_p, D01summ_coef_p, D03summ_coef_p, D08summ_coef_p, ncol = 2, nrow = 2) 
ggarrange(D01acc_plot, D01acc_plot, D03acc_plot, D08acc_plot, ncol = 2, nrow = 2) 
