library(ggplot2)
# Numerical illustration --------------------------------------------------

plot1_po_pn <- plot1_po_n <- plot1_o_pn <- plot1_o_n <- NULL
set.seed(30) # 4  ou 11
g_deviance <- function(nu, data_in) {
  logLik(glm(y~.,family = binomial(Gosset(nu)),data = data_in, maxit = 100))
}

num_exp <- function(n_obs,s,n_out,df_stu,seeed,out_y="meth2"){
  set.seed(seeed)
  
  n_obs <- n_obs*1.3
  
  alpha <- 0
  beta1 <- c(1, rep(0, 1))
  Sigma <- matrix(0, nrow = 2, ncol = 2)
  diag(Sigma) <- 1
  Ex <- cbind(MASS::mvrnorm(n_obs, rep(0, 2), Sigma))
  eta1 <- alpha + as.matrix(Ex) %*% (beta1)
  h_x_beta <- S(eta1,s)
  y <- rbinom(n = n_obs,size = 1,prob = h_x_beta)  
  dat2 <- data.frame(y,Ex)
  dat2$y <- as.factor(dat2$y)
  dat3 <- dat2
  
  n_out_f <- n_out
  Sigma2 <- matrix(0, nrow = 2, ncol = 2)
  diag(Sigma2) <- 0.5
  if(n_out_f > 0){
    if(n_out_f==1){
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = 2, rep(2, 2), Sigma2)),y=rep(0,n_out_f)) 
      dat3 <- rbind(dat2,outliers[1,])
    }else{
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = n_out_f, rep(2, 2), Sigma2)),y=rep(0,n_out_f)) 
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
      log_nu[i] <- logLik(glm(formula = y ~. , data = dat3, family = binomial(link = Gosset(nu)),maxit = 100))
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
A005acc_stu <- A005acc_log <- A005df_stu <- A005coefs_b <- A005coefs_f <- A005coefs <- A005AICsteplog <- A005BICstepstu <- A005AICstepstu <-  A005BICsteplog <- A005max_grille <- A005df_stu_opt <- NULL
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
    
    mtry <- try(num_exp(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable)
    A005coefs <- rbind(A005coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A005acc_log <- rbind(A005acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    A005acc_stu <- rbind(A005acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    A005AICsteplog <- rbind(A005AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A005AICstepstu <- rbind(A005AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A005BICsteplog <- rbind(A005BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A005BICstepstu <- rbind(A005BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # A005df_stu_opt <- rbind(A005df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    A005df_stu <- rbind(A005df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(A005coefs)
dat5 <- A005coefs[,-ncol(A005coefs)]
dat5[dat5==0] <- NA
A005coefs[,-ncol(A005coefs)] <- dat5
colnames(A005coefs)[ncol(A005coefs)] <- "out"


save(A005coefs, file = "data/A005coefs.RData")
# load("$data/A005coefs.RData")

A005summ_coef <- A005coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A005summ_coef)


A005summ_coef_p <- ggplot(data = A005summ_coef, aes(x=out,color=link)) +
  geom_line(aes(y=X2, linetype = "dashed")) +
  geom_line(aes(y=X1,linetype = "solid")) +
  geom_point(aes(y=X1))+
  geom_point(aes(y=X2))+
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.05$")

A005summ_coef_p

tikz(file = "tex/A005summ_coef_p2.tex", standAlone=F, width = 7.5, height = 4.5)
# endoffile <- dev.off()


A005acc_logit <- as.data.frame(A005acc_log)
colnames(A005acc_logit)[1] <- "loglik"
A005acc_logit$link <- "logit"

A005acc_stu <- as.data.frame(A005acc_stu)
colnames(A005acc_stu)[1] <- "loglik"
A005acc_stu$link <- "stu"

A005accliks <- rbind(A005acc_logit,A005acc_stu)
head(A005accliks)
colnames(A005accliks) <- c("loglik","outlier","link")

save(A005accliks, file = "data/A005accliks.RData")

A005acc_plot <- ggplot(data = A005accliks[A005accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.05")

A005acc_plot

# 0.1 --------------------------------------------------------------------

nsim <- 100
A01acc_stu <- A01acc_log <- A01df_stu <- A01coefs_b <- A01coefs_f <- A01coefs <- A01AICsteplog <- A01BICstepstu <- A01AICstepstu <-  A01BICsteplog <- A01max_grille <- A01df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable)
    A01coefs <- rbind(A01coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A01acc_log <- rbind(A01acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    A01acc_stu <- rbind(A01acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    A01AICsteplog <- rbind(A01AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A01AICstepstu <- rbind(A01AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A01BICsteplog <- rbind(A01BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A01BICstepstu <- rbind(A01BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # A01df_stu_opt <- rbind(A01df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    A01df_stu <- rbind(A01df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(A01coefs)
dat5 <- A01coefs[,-ncol(A01coefs)]
dat5[dat5==0] <- NA
A01coefs[,-ncol(A01coefs)] <- dat5
colnames(A01coefs)[ncol(A01coefs)] <- "out"

save(A01coefs, file = "data/A01coefs.RData")

A01summ_coef <- A01coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A01summ_coef)

A01summ_coef_p <- ggplot(data = A01summ_coef, aes(x=out,color=link)) +
  geom_line(aes(y=X2, linetype = "dashed")) +
  geom_line(aes(y=X1,linetype = "solid")) +
  geom_point(aes(y=X1))+
  geom_point(aes(y=X2))+
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.1$")

A01summ_coef_p

A01acc_logit <- as.data.frame(A01acc_log)
colnames(A01acc_logit)[1] <- "loglik"
A01acc_logit$link <- "logit"

A01acc_stu <- as.data.frame(A01acc_stu)
colnames(A01acc_stu)[1] <- "loglik"
A01acc_stu$link <- "stu"

A01accliks <- rbind(A01acc_logit,A01acc_stu)
head(A01accliks)
colnames(A01accliks) <- c("loglik","outlier","link")

save(A01accliks, file = "data/A01accliks.RData")

A01acc_plot <- ggplot(data = A01accliks[A01accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.1")

A01acc_plot

# 0.3 --------------------------------------------------------------------

nsim <- 100
A03acc_stu <- A03acc_log <- A03df_stu <- A03coefs_b <- A03coefs_f <- A03coefs <- A03AICsteplog <- A03BICstepstu <- A03AICstepstu <-  A03BICsteplog <- A03max_grille <- A03df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable)
    A03coefs <- rbind(A03coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A03acc_log <- rbind(A03acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    A03acc_stu <- rbind(A03acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    A03AICsteplog <- rbind(A03AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A03AICstepstu <- rbind(A03AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A03BICsteplog <- rbind(A03BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A03BICstepstu <- rbind(A03BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # A03df_stu_opt <- rbind(A03df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    A03df_stu <- rbind(A03df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(A03coefs)
dat5 <- A03coefs[,-ncol(A03coefs)]
dat5[dat5==0] <- NA
A03coefs[,-ncol(A03coefs)] <- dat5
colnames(A03coefs)[ncol(A03coefs)] <- "out"

save(A03coefs, file = "data/A03coefs.RData")

A03summ_coef <- A03coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A03summ_coef)

A03summ_coef_p <- ggplot(data = A03summ_coef, aes(x=out,color=link)) +
  geom_line(aes(y=X2, linetype = "dashed")) +
  geom_line(aes(y=X1,linetype = "solid")) +
  geom_point(aes(y=X1))+
  geom_point(aes(y=X2))+
  xlab(" of outliers") + ylab(" of var. selection") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.3$")

A03summ_coef_p

A03acc_logit <- as.data.frame(A03acc_log)
colnames(A03acc_logit)[1] <- "loglik"
A03acc_logit$link <- "logit"

A03acc_stu <- as.data.frame(A03acc_stu)
colnames(A03acc_stu)[1] <- "loglik"
A03acc_stu$link <- "stu"

A03accliks <- rbind(A03acc_logit,A03acc_stu)
head(A03accliks)
colnames(A03accliks) <- c("loglik","outlier","link")

save(A03accliks, file = "data/A03accliks.RData")

A03acc_plot <- ggplot(data = A03accliks[A03accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.3")

A03acc_plot

# 0.8 --------------------------------------------------------------------

nsim <- 100
A08acc_stu <- A08acc_log <- A08df_stu <- A08coefs_b <- A08coefs_f <- A08coefs <- A08AICsteplog <- A08BICstepstu <- A08AICstepstu <-  A08BICsteplog <- A08max_grille <- A08df_stu_opt <- NULL
variable <- 1
variable_e <- 1
for (outlier in 0:20) {
  # outlier <- 0
  for (variable_e in 1:nsim) {
    variable_e <- variable_e+1
    variable <- variable+1
    
    print(paste0(outlier,"A",variable))
    
    mtry <- try(num_exp(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable),
      silent=TRUE)
    } 
    mode1 <- mtry
    # mode1 <- num_exp(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable)
    A08coefs <- rbind(A08coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    A08acc_log <- rbind(A08acc_log,cbind(mode1[["acc_log"]],rbind(outlier)))
    A08acc_stu <- rbind(A08acc_stu,cbind(mode1[["acc_stu"]],rbind(outlier)))
    A08AICsteplog <- rbind(A08AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A08AICstepstu <- rbind(A08AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A08BICsteplog <- rbind(A08BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A08BICstepstu <- rbind(A08BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    # max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    # A08df_stu_opt <- rbind(A08df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    A08df_stu <- rbind(A08df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

head(A08coefs)
dat5 <- A08coefs[,-ncol(A08coefs)]
dat5[dat5==0] <- NA
A08coefs[,-ncol(A08coefs)] <- dat5
colnames(A08coefs)[ncol(A08coefs)] <- "out"

save(A08coefs, file = "data/A08coefs.RData")

A08summ_coef <- A08coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A08summ_coef)

A08summ_coef_p <- ggplot(data = A08summ_coef, aes(x=out,color=link)) +
  geom_line(aes(y=X2, linetype = "dashed")) +
  geom_line(aes(y=X1,linetype = "solid")) +
  geom_point(aes(y=X1))+
  geom_point(aes(y=X2))+
  xlab(" of outliers") + ylab("of var. selection") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.8$")

A08summ_coef_p

A08acc_logit <- as.data.frame(A08acc_log)
colnames(A08acc_logit)[1] <- "loglik"
A08acc_logit$link <- "logit"

A08acc_stu <- as.data.frame(A08acc_stu)
colnames(A08acc_stu)[1] <- "loglik"
A08acc_stu$link <- "stu"

A08accliks <- rbind(A08acc_logit,A08acc_stu)
head(A08accliks)
colnames(A08accliks) <- c("loglik","outlier","link")

save(A08accliks, file = "data/A08accliks.RData")

A08acc_plot <- ggplot(data = A08accliks[A08accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("") + ylab("") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("0.8")

A08acc_plot

# tikz(file = "tex/Av_coef_2_Noisy.tex", standAlone=F, width = 7.5, height = 4.5)
# endoffile <- dev.off()

ggarrange(A005summ_coef_p, A01summ_coef_p, A03summ_coef_p, A08summ_coef_p, ncol = 2, nrow = 2) 
ggarrange(A005acc_plot, A01acc_plot, A03acc_plot, A08acc_plot, ncol = 2, nrow = 2) 
