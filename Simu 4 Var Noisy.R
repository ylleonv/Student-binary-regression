library(ggplot2)
# Numerical illustration --------------------------------------------------

plot1_po_pn <- plot1_po_n <- plot1_o_pn <- plot1_o_n <- NULL
set.seed(30) # 4  ou 11
g_deviance <- function(nu, data_in) {
  logLik(glm(y~.,family = binomial(Gosset(nu)),data = data_in, maxit = 100))
}

num_exp4 <- function(n_obs,s,n_out,df_stu,seeed,out_y="meth2"){
  set.seed(seeed)
  
  n_obs <- n_obs*1.3
  
  alpha <- 0
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
  
  n_out_f <- n_out
  Sigma2 <- matrix(0, nrow = 4, ncol = 4)
  diag(Sigma2) <- 0.5
  if(n_out_f > 0){
    if(n_out_f==1){
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = 2, c(1.6,-1.2,2,2), Sigma2)),y=rep(0,n_out_f)) 
      dat3 <- rbind(dat2,outliers[1,])
    }else{
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = n_out_f, c(1.6,-1.2,2,2), Sigma2)),y=rep(0,n_out_f)) 
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
B005acc_stu <- B005acc_log <- B005df_stu <- B005coefs_b <- B005coefs_f <- B005coefs <- B005AICsteplog <- B005BICstepstu <- B005AICstepstu <-  B005BICsteplog <- B005max_grille <- B005df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp4(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp4(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp4(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable)
    B005coefs <- rbind(B005coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    B005acc_log <- rbind(B005acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    B005acc_stu <- rbind(B005acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    B005AICsteplog <- rbind(B005AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    B005AICstepstu <- rbind(B005AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    B005BICsteplog <- rbind(B005BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    B005BICstepstu <- rbind(B005BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # B005df_stu_opt <- rbind(B005df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    B005df_stu <- rbind(B005df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(B005coefs)
dat5 <- B005coefs[,-ncol(B005coefs)]
dat5[dat5==0] <- NA
B005coefs[,-ncol(B005coefs)] <- dat5
colnames(B005coefs)[ncol(B005coefs)] <- "out"

B005summ_coef <- B005coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(B005summ_coef)

save(B005summ_coef, file = "data/Simu 4 Var Noisy/B005summ_coef.RData")

B005_SUM <- cbind(B005summ_coef[,1:2],
  exva = rowMeans(B005summ_coef[,4:5]),
  novar = rowMeans(B005summ_coef[,6:7]))

B005summ_coef_p <- ggplot(data = B005_SUM, aes(x=out,color=link)) +
  geom_line(aes(y=novar, linetype = "dashed")) +
  geom_line(aes(y=exva,linetype = "solid")) +
  geom_point(aes(y=novar))+
  geom_point(aes(y=exva))+
  xlab("of outliers") + ylab("of var. selection") +
  scale_color_manual("", aesthetics = "colour",labels = c("Logistic", "Student"), values = c("#E69F00", "#00AFBB"))+
  scale_linetype_manual("", labels = c("Discriminant", "Noisy"), values = c("dashed", "solid"))+
  theme(
    strip.background = element_rect(
      fill="yellow"
    ),
    panel.background = element_rect(fill = "yellow", colour = "grey50"),
    plot.background = element_rect(fill = "yellow", colour = "yellow"),
    legend.background = element_rect(fill = "yellow"),
    legend.key = element_rect(fill = NA),
    legend.position = "top"
  )  + ggtitle("$d=0.05$")

B005summ_coef_p

tikz(file = "tex/Label_variable.tex", standAlone=F, width = 7.5, height = 4.5)
endoffile <- dev.off()

B005acc_logit <- as.data.frame(B005acc_log)
colnames(B005acc_logit)[1] <- "loglik"
B005acc_logit$link <- "logit"

B005acc_stu <- as.data.frame(B005acc_stu)
colnames(B005acc_stu)[1] <- "loglik"
B005acc_stu$link <- "stu"

B005accliks <- rbind(B005acc_logit,B005acc_stu)
head(B005accliks)
colnames(B005accliks) <- c("loglik","outlier","link")

save(B005accliks, file = "data/Simu 4 Var Noisy/B005accliks.RData")

B005acc_plot <- ggplot(data = B005accliks[B005accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("of var. selection") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.05")

B005acc_plot

# 0.1 --------------------------------------------------------------------

nsim <- 100
B01acc_stu <- B01acc_log <- B01df_stu <- B01coefs_b <- B01coefs_f <- B01coefs <- B01AICsteplog <- B01BICstepstu <- B01AICstepstu <-  B01BICsteplog <- B01max_grille <- B01df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp4(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp4(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp4(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable)
    B01coefs <- rbind(B01coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    B01acc_log <- rbind(B01acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    B01acc_stu <- rbind(B01acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    B01AICsteplog <- rbind(B01AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    B01AICstepstu <- rbind(B01AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    B01BICsteplog <- rbind(B01BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    B01BICstepstu <- rbind(B01BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # B01df_stu_opt <- rbind(B01df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    B01df_stu <- rbind(B01df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(B01coefs)
dat5 <- B01coefs[,-ncol(B01coefs)]
dat5[dat5==0] <- NA
B01coefs[,-ncol(B01coefs)] <- dat5
colnames(B01coefs)[ncol(B01coefs)] <- "out"

B01summ_coef <- B01coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(B01summ_coef)

save(B01summ_coef, file = "data/Simu 4 Var Noisy/B01summ_coef.RData")

B01_SUM <- cbind(B01summ_coef[,1:2],
  exva = rowMeans(B01summ_coef[,4:5]),
  novar = rowMeans(B01summ_coef[,6:7]))

B01summ_coef_p <- ggplot(data = B01_SUM, aes(x=out,color=link)) +
  geom_line(aes(y=novar, linetype = "dashed")) +
  geom_line(aes(y=exva,linetype = "solid")) +
  geom_point(aes(y=novar))+
  geom_point(aes(y=exva))+
  xlab("of outliers") + ylab("of var. selection") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.1$")

B01summ_coef_p

B01acc_logit <- as.data.frame(B01acc_log)
colnames(B01acc_logit)[1] <- "loglik"
B01acc_logit$link <- "logit"

B01acc_stu <- as.data.frame(B01acc_stu)
colnames(B01acc_stu)[1] <- "loglik"
B01acc_stu$link <- "stu"

B01accliks <- rbind(B01acc_logit,B01acc_stu)
head(B01accliks)
colnames(B01accliks) <- c("loglik","outlier","link")

save(B01accliks, file = "data/Simu 4 Var Noisy/B01accliks.RData")

B01acc_plot <- ggplot(data = B01accliks[B01accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("of var. selection") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.1")

B01acc_plot

# 0.3 --------------------------------------------------------------------

nsim <- 100
B03acc_stu <- B03acc_log <- B03df_stu <- B03coefs_b <- B03coefs_f <- B03coefs <- B03AICsteplog <- B03BICstepstu <- B03AICstepstu <-  B03BICsteplog <- B03max_grille <- B03df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp4(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp4(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp4(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable)
    B03coefs <- rbind(B03coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    B03acc_log <- rbind(B03acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    B03acc_stu <- rbind(B03acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    B03AICsteplog <- rbind(B03AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    B03AICstepstu <- rbind(B03AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    B03BICsteplog <- rbind(B03BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    B03BICstepstu <- rbind(B03BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # B03df_stu_opt <- rbind(B03df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    B03df_stu <- rbind(B03df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(B03coefs)
dat5 <- B03coefs[,-ncol(B03coefs)]
dat5[dat5==0] <- NA
B03coefs[,-ncol(B03coefs)] <- dat5
colnames(B03coefs)[ncol(B03coefs)] <- "out"

B03summ_coef <- B03coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(B03summ_coef)

save(B03summ_coef, file = "data/Simu 4 Var Noisy/B03summ_coef.RData")

B03_SUM <- cbind(B03summ_coef[,1:2],
  exva = rowMeans(B03summ_coef[,4:5]),
  novar = rowMeans(B03summ_coef[,6:7]))

B03summ_coef_p <- ggplot(data = B03_SUM, aes(x=out,color=link)) +
  geom_line(aes(y=novar, linetype = "dashed")) +
  geom_line(aes(y=exva,linetype = "solid")) +
  geom_point(aes(y=novar))+
  geom_point(aes(y=exva))+
  xlab("of outliers") + ylab("of var. selection") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.3$")

B03summ_coef_p

B03acc_logit <- as.data.frame(B03acc_log)
colnames(B03acc_logit)[1] <- "loglik"
B03acc_logit$link <- "logit"

B03acc_stu <- as.data.frame(B03acc_stu)
colnames(B03acc_stu)[1] <- "loglik"
B03acc_stu$link <- "stu"

B03accliks <- rbind(B03acc_logit,B03acc_stu)
head(B03accliks)
colnames(B03accliks) <- c("loglik","outlier","link")

save(B03accliks, file = "data/Simu 4 Var Noisy/B03accliks.RData")

B03acc_plot <- ggplot(data = B03accliks[B03accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("of var. selection") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.3")

B03acc_plot

# 0.8 --------------------------------------------------------------------

nsim <- 100
B08acc_stu <- B08acc_log <- B08df_stu <- B08coefs_b <- B08coefs_f <- B08coefs <- B08AICsteplog <- B08BICstepstu <- B08AICstepstu <-  B08BICsteplog <- B08max_grille <- B08df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp4(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp4(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp4(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable)
    B08coefs <- rbind(B08coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    B08acc_log <- rbind(B08acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    B08acc_stu <- rbind(B08acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    B08AICsteplog <- rbind(B08AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    B08AICstepstu <- rbind(B08AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    B08BICsteplog <- rbind(B08BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    B08BICstepstu <- rbind(B08BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # B08df_stu_opt <- rbind(B08df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    B08df_stu <- rbind(B08df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(B08coefs)
dat5 <- B08coefs[,-ncol(B08coefs)]
dat5[dat5==0] <- NA
B08coefs[,-ncol(B08coefs)] <- dat5
colnames(B08coefs)[ncol(B08coefs)] <- "out"

B08summ_coef <- B08coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(B08summ_coef)

save(B08summ_coef, file = "data/Simu 4 Var Noisy/B08summ_coef.RData")

B08_SUM <- cbind(B08summ_coef[,1:2],
  exva = rowMeans(B08summ_coef[,4:5]),
  novar = rowMeans(B08summ_coef[,6:7]))

B08summ_coef_p <- ggplot(data = B08_SUM, aes(x=out,color=link)) +
  geom_line(aes(y=novar, linetype = "dashed")) +
  geom_line(aes(y=exva,linetype = "solid")) +
  geom_point(aes(y=novar))+
  geom_point(aes(y=exva))+
  xlab("of outliers") + ylab("of var. selection") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.8$")

B08summ_coef_p

B08acc_logit <- as.data.frame(B08acc_log)
colnames(B08acc_logit)[1] <- "loglik"
B08acc_logit$link <- "logit"

B08acc_stu <- as.data.frame(B08acc_stu)
colnames(B08acc_stu)[1] <- "loglik"
B08acc_stu$link <- "stu"

B08accliks <- rbind(B08acc_logit,B08acc_stu)
head(B08accliks)
colnames(B08accliks) <- c("loglik","outlier","link")

save(B08accliks, file = "data/Simu 4 Var Noisy/B08accliks.RData")

B08acc_plot <- ggplot(data = B08accliks[B08accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("of var. selection") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.8")

B08acc_plot

# tikz(file = "tex/Av_coef_4_Noisy.tex", standAlone=F, width = 7.5, height = 4.5)
# endoffile <- dev.off()

ggarrange(B005summ_coef_p, B01summ_coef_p, B03summ_coef_p, B08summ_coef_p, ncol = 2, nrow = 2) 
ggarrange(B005acc_plot, B01acc_plot, B03acc_plot, B08acc_plot, ncol = 2, nrow = 2) 
