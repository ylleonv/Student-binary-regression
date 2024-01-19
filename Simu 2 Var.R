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
  # df_stu <- 0.8
  mod_log <- glm(formula = y ~ X1 , data = dat3, family = binomial(link = "logit"))
  mod_stu <- glm(formula = y ~ X1 , data = dat3, family = binomial(link = Gosset(df_stu)))
  plot1_po_pn <- ggplot(dat3,aes(x=X1,y=y,color=y))+
    geom_point()+
    scale_color_manual("", values = c("#E69F00", "#00AFBB"))+
    geom_vline(xintercept = -coef(mod_log)[1]/coef(mod_log)[2])+
    geom_vline(xintercept = -coef(mod_stu)[1]/coef(mod_stu)[2], linetype = "dashed")+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
      legend.key = element_rect(fill = "white"),
      legend.position = "none")+
    ylab("")+xlab("")
  plot1_po_pn
  mod_log <- glm(formula = y ~ X1 + X2, data = dat3, family = binomial(link = "logit"))
  mod_stu <- glm(formula = y ~ X1 + X2, data = dat3, family = binomial(link = Gosset(df_stu)))
  
  plot1_po_n <- ggplot(dat3,aes(x=X1,y=X2,color=y))+
    geom_point()+
    scale_color_manual("", values = c("#E69F00", "#00AFBB"))+
    geom_vline(xintercept = -coef(mod_log)[1]/coef(mod_log)[2])+
    # ylim(-3,3)+
    geom_vline(xintercept = -coef(mod_stu)[1]/coef(mod_stu)[2], linetype = "dashed")+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
      legend.key = element_rect(fill = "white"),
      legend.position = "none")+
    ylab("")+xlab("")
  plot1_po_n
  
  n_out_f <- n_out
  Sigma2 <- matrix(0, nrow = 2, ncol = 2)
  diag(Sigma2) <- 0.25
  if(n_out_f > 0){
    if(n_out_f==1){
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = 2, rep(2, 2), Sigma2)),y=rep(0,n_out_f)) 
      dat3 <- rbind(dat2,outliers[1,])
    }else{
      outliers <- data.frame(data.frame(MASS::mvrnorm(n = n_out_f, rep(2, 2), Sigma2)),y=rep(0,n_out_f)) 
      dat3 <- rbind(dat2,outliers)
    }
  }
  
  nu_expl <- c(seq(0.3,1,0.1),2:8)
  log_nu <- NULL
  i <- 1
  for (nu in c(seq(0.3,1,0.1),2:8)) {
    log_nu[i] <- logLik(glm(formula = y ~ X1 + X2, data = dat3, family = binomial(link = Gosset(nu))))
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
  
  
  
  mod_log <- glm(formula = y ~ X1 , data = dat3, family = binomial(link = "logit"))
  mod_stu <- glm(formula = y ~ X1, data = dat3, family = binomial(link = Gosset(df_stu)))
  plot1_o_pn <- ggplot(dat3,aes(x=X1,y=y,color=y))+
    geom_point()+
    scale_color_manual("", values = c("#E69F00", "#00AFBB"))+
    geom_vline(xintercept = -coef(mod_log)[1]/coef(mod_log)[2])+
    geom_vline(xintercept = -coef(mod_stu)[1]/coef(mod_stu)[2], linetype = "dashed")+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
      legend.key = element_rect(fill = "white"),
      legend.position = "none")+
    ylab("")+xlab("")
  
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
  
  plot1_o_n <- ggplot(dat3,aes(x=X1,y=X2,color=y))+
    geom_point()+
    scale_color_manual("", values = c("#E69F00", "#00AFBB"))+
    stat_function(fun = function(x) (-coef(mod_log_n1)[1]/coef(mod_log_n1)[3])-((coef(mod_log_n1)[2]/coef(mod_log_n1)[3])*x),
      linetype = "solid", color = "black")+
    stat_function(fun = function(x) (-coef(mod_stu_n1)[1]/coef(mod_stu_n1)[3])-((coef(mod_stu_n1)[2]/coef(mod_stu_n1)[3])*x),
      linetype = "dashed", color = "black")+
    ylim(-3,3)+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
      legend.key = element_rect(fill = "white"),
      legend.position = "none")
  
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
    # coef_sel_f = dat_b_f,
    coef_sel_b = dat_b_b,
    max_grille = max_grille,df_stu_opt = df_stu_opt,
    df_stu = df_stu
    ))
}



# 0.05 --------------------------------------------------------------------

nsim <- 100
df_stu <- coefs_b <- coefs_f <- coefs <- AICsteplog <- BICstepstu <- AICstepstu <-  BICsteplog <- max_grille <- df_stu_opt <- NULL
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp(n_obs = 100,s = 0.05,n_out = outlier,df_stu = 0.8,seeed = variable)
    coefs <- rbind(coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    # coefs_f <- rbind(coefs_f,cbind(mode1[["coef_sel_f"]],rbind(outlier,outlier)))
    coefs_b <- rbind(coefs_b,cbind(mode1[["coef_sel_b"]],rbind(outlier,outlier)))
    AICsteplog <- rbind(AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    AICstepstu <- rbind(AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    BICsteplog <- rbind(BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    BICstepstu <- rbind(BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    max_grille <- rbind(max_grille,cbind(mode1[["max_grille"]],outlier))
    df_stu_opt <- rbind(df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    df_stu <- rbind(df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

A05coefs_f <- coefs_f
A05coefs_b <- coefs_b

degree_info <- data.frame(cbind(max_grille, df_stu_opt, df_stu[-nrow(df_stu)]))
colnames(degree_info) <- c("maxgri","outlier","maxalgo","out2","selected")
degree_info$sel_algo <- c(degree_info$maxalgo == degree_info$selected)
head(degree_info)

plot_zoom <- ggplot(data = degree_info, mapping = aes(x=maxgri,y=selected, shape = sel_algo, colour = outlier))+
  # geom_point()+
  coord_cartesian(xlim=c(0.25,1), ylim = c(0.2,1))+
  geom_jitter()+
  xlab("v obtained from the grid") + ylab("v obtained with optimization algo") +
  scale_shape_manual("",labels = c("v = 8", "Selected"), values = c(1, 4))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    # legend.position = "none"
  )  + ggtitle("")

plot_com <- ggplot(data = degree_info, mapping = aes(x=maxgri,y=selected, shape = sel_algo, colour = outlier))+
  # geom_point()+
  # coord_cartesian(xlim=c(0.25,1), ylim = c(0.2,1))+
  geom_jitter()+
  xlab("v obtained from the grid") + ylab("v obtained with optimization algo") +
  scale_shape_manual("",labels = c("$v = 8$", "$Selected$"), values = c(1, 4))+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("")

ggarrange(plot_com, plot_zoom, ncol = 2, nrow = 1)

ggplot(data = degree_info, mapping = aes(x=factor(outlier),y=selected,colour=sel_algo))+
  geom_jitter(height = 0.2)
  # coord_cartesian(xlim=c(0.25,1))+
  # geom_boxplot()


dat5 <- coefs[,-5]
dat5[dat5==0] <- NA
coefs[,-5] <- dat5
colnames(coefs)[5] <- "out"
summ_coef <- coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(summ_coef)

BIC_PLOT <- ggplot(data = summ_coef, aes(x=out,color=link)) +
  geom_line(aes(y=X2, linetype = "dashed")) +
  geom_line(aes(y=X1,linetype = "solid")) +
  geom_point(aes(y=X1))+
  geom_point(aes(y=X2))+
  xlab("") + ylab("") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.05")

BIC_PLOT

# 0.1 --------------------------------------------------------------------

nsim <- 100
A1df_stu <- A1coefs_b <- A1coefs_f <- A1coefs <- A1AICsteplog <- A1BICstepstu <- A1AICstepstu <-  A1BICsteplog <- A1max_grille <- A1df_stu_opt <- NULL
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp(n_obs = 100,s = 0.1,n_out = outlier,df_stu = 0.8,seeed = variable)
    A1coefs <- rbind(A1coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    # A1coefs_f <- rbind(A1coefs_f,cbind(mode1[["coef_sel_f"]],rbind(outlier,outlier)))
    A1coefs_b <- rbind(A1coefs_b,cbind(mode1[["coef_sel_b"]],rbind(outlier,outlier)))
    A1AICsteplog <- rbind(A1AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A1AICstepstu <- rbind(A1AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A1BICsteplog <- rbind(A1BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A1BICstepstu <- rbind(A1BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    A1max_grille <- rbind(A1max_grille,cbind(mode1[["max_grille"]],outlier))
    A1df_stu_opt <- rbind(A1df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    A1df_stu <- rbind(A1df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

cbind(A1max_grille, A1df_stu_opt)

A1dat5 <- A1coefs[,-5]
A1dat5[A1dat5==0] <- NA
A1coefs[,-5] <- A1dat5
colnames(A1coefs)[5] <- "out"
A1summ_coef <- A1coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A1summ_coef)

A1BIC_PLOT <- ggplot(data = A1summ_coef, aes(x=out,color=link)) +
  geom_line(aes(y=X2, linetype = "dashed")) +
  geom_line(aes(y=X1,linetype = "solid")) +
  geom_point(aes(y=X1))+
  geom_point(aes(y=X2))+
  xlab("") + ylab("") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.1")

A1BIC_PLOT

A1dat5_coef_b <- A1coefs_b
A1dat5_coef_b[A1dat5_coef_b==0] <- NA
colnames(A1dat5_coef_b)[ncol(A1dat5_coef_b)] <- "out"
A1dat5_coef_b <- A1dat5_coef_b %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A1dat5_coef_b)

A1BIC_PLOT_B <- ggplot(data = A1dat5_coef_b, aes(x=out,color=link)) +
  geom_line(aes(y=X2, linetype = "dashed")) +
  geom_line(aes(y=X1,linetype = "solid")) +
  geom_point(aes(y=X1))+
  geom_point(aes(y=X2))+
  xlab("") + ylab("") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.1")

A1BIC_PLOT_B

A1dat5_coef_F <- A1coefs_f
A1dat5_coef_F[A1dat5_coef_F==0] <- NA
colnames(A1dat5_coef_F)[ncol(A1dat5_coef_F)] <- "out"
A1dat5_coef_F <- A1dat5_coef_F %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A1dat5_coef_F)

A1BIC_PLOT_F <- ggplot(data = A1dat5_coef_F, aes(x=out,color=link)) +
  geom_line(aes(y=X2, linetype = "dashed")) +
  geom_line(aes(y=X1,linetype = "solid")) +
  geom_point(aes(y=X1))+
  geom_point(aes(y=X2))+
  xlab("") + ylab("") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.1")

A1BIC_PLOT_F

# 0.3 --------------------------------------------------------------------

nsim <- 100
A3df_stu <- A3coefs_b <- A3coefs_f <- A3coefs <- A3AICsteplog <- A3BICstepstu <- A3AICstepstu <-  A3BICsteplog <- A3max_grille <- A3df_stu_opt <- NULL
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp(n_obs = 100,s = 0.3,n_out = outlier,df_stu = 0.8,seeed = variable)
    A3coefs <- rbind(A3coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    # A3coefs_f <- rbind(A3coefs_f,cbind(mode1[["coef_sel_f"]],rbind(outlier,outlier)))
    A3coefs_b <- rbind(A3coefs_b,cbind(mode1[["coef_sel_b"]],rbind(outlier,outlier)))
    A3AICsteplog <- rbind(A3AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A3AICstepstu <- rbind(A3AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A3BICsteplog <- rbind(A3BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A3BICstepstu <- rbind(A3BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    A3max_grille <- rbind(A3max_grille,cbind(mode1[["max_grille"]],outlier))
    A3df_stu_opt <- rbind(A3df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    A3df_stu <- rbind(A3df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

cbind(A3max_grille, A3df_stu_opt)

A3dat5 <- A3coefs[,-5]
A3dat5[A3dat5==0] <- NA
A3coefs[,-5] <- A3dat5
colnames(A3coefs)[5] <- "out"
A3summ_coef <- A3coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A3summ_coef)

A3BIC_PLOT <- ggplot(data = A3summ_coef, aes(x=out,color=link)) +
  geom_line(aes(y=X2, linetype = "dashed")) +
  geom_line(aes(y=X1,linetype = "solid")) +
  geom_point(aes(y=X1))+
  geom_point(aes(y=X2))+
  xlab("") + ylab("") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.3")

A3BIC_PLOT

# 0.8 --------------------------------------------------------------------

nsim <- 100
A8df_stu <- A8coefs_b <- A8coefs_f <- A8coefs <- A8AICsteplog <- A8BICstepstu <- A8AICstepstu <-  A8BICsteplog <- A8max_grille <- A8df_stu_opt <- NULL
for (outlier in 1:20) {
  # outlier <- 0
  for (variable in 1:nsim) {
    # variable <- 1
    mode1 <- num_exp(n_obs = 100,s = 0.8,n_out = outlier,df_stu = 0.8,seeed = variable)
    A8coefs <- rbind(A8coefs,cbind(mode1[["coef_sel"]],rbind(outlier,outlier)))
    # A8coefs_f <- rbind(A8coefs_f,cbind(mode1[["coef_sel_f"]],rbind(outlier,outlier)))
    A8coefs_b <- rbind(A8coefs_b,cbind(mode1[["coef_sel_b"]],rbind(outlier,outlier)))
    A8AICsteplog <- rbind(A8AICsteplog,cbind(mode1[["AIC_step_log"]],outlier))
    A8AICstepstu <- rbind(A8AICstepstu,cbind(mode1[["AIC_step_stu"]],outlier))
    A8BICsteplog <- rbind(A8BICsteplog,cbind(mode1[["BIC_step_log"]],outlier))
    A8BICstepstu <- rbind(A8BICstepstu,cbind(mode1[["BIC_step_stu"]],outlier))
    A8max_grille <- rbind(A8max_grille,cbind(mode1[["max_grille"]],outlier))
    A8df_stu_opt <- rbind(A8df_stu_opt,cbind(mode1[["df_stu_opt"]],outlier))
    A8df_stu <- rbind(A8df_stu,cbind(mode1[["df_stu"]], outlier))
  }
  print(outlier)
}

cbind(A8max_grille, A8df_stu_opt)

A8dat5 <- A8coefs[,-5]
A8dat5[A8dat5==0] <- NA
A8coefs[,-5] <- A8dat5
colnames(A8coefs)[5] <- "out"
A8summ_coef <- A8coefs %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A8summ_coef)

A8BIC_PLOT <- ggplot(data = A8summ_coef, aes(x=out,color=link)) +
  geom_line(aes(y=X2, linetype = "dashed")) +
  geom_line(aes(y=X1,linetype = "solid")) +
  geom_point(aes(y=X1))+
  geom_point(aes(y=X2))+
  xlab("") + ylab("") +
  scale_color_manual("", aesthetics = "colour",labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("s=0.8")

A8BIC_PLOT


ggarrange(BIC_PLOT, A1BIC_PLOT, A3BIC_PLOT, A8BIC_PLOT, ncol = 2, nrow = 2) 


# Graficos ----------------------------------------------------------------

A05dat5 <- A05coefs_b
A05dat5[A05dat5==0] <- NA
A05coefs_b <- A05dat5
colnames(A05coefs_b)[ncol(A05coefs_b)] <- "out"
A05summ_coef <- A05coefs_b %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A05summ_coef)

A05_SUM <- cbind(A05summ_coef[,1:2],
  exva = rowMeans(A05summ_coef[,4]),
  novar = rowMeans(A05summ_coef[,5]))

A05BIC_PLOT <- ggplot(data = A05_SUM, aes(x=out,color=link)) +
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

A05BIC_PLOT

A1dat5 <- A1coefs_b
A1dat5[A1dat5==0] <- NA
A1coefs_b <- A1dat5
colnames(A1coefs_b)[ncol(A1coefs_b)] <- "out"
A1summ_coef <- A1coefs_b %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A1summ_coef)

A1_SUM <- cbind(A1summ_coef[,1:2],
  exva = rowMeans(A1summ_coef[,4]),
  novar = rowMeans(A1summ_coef[,5]))

A1BIC_PLOT <- ggplot(data = A1_SUM, aes(x=out,color=link)) +
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

A1BIC_PLOT

A3dat5 <- A3coefs_b
A3dat5[A3dat5==0] <- NA
A3coefs_b <- A3dat5
colnames(A3coefs_b)[ncol(A3coefs_b)] <- "out"
A3summ_coef <- A3coefs_b %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A3summ_coef)

A3_SUM <- cbind(A3summ_coef[,1:2],
  exva = rowMeans(A3summ_coef[,4]),
  novar = rowMeans(A3summ_coef[,5]))

A3BIC_PLOT <- ggplot(data = A3_SUM, aes(x=out,color=link)) +
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

A3BIC_PLOT

A8dat5 <- A8coefs_b
A8dat5[A8dat5==0] <- NA
A8coefs_b <- A8dat5
colnames(A8coefs_b)[ncol(A8coefs_b)] <- "out"
A8summ_coef <- A8coefs_b %>% 
  group_by(out,link) %>% 
  dplyr::summarise_all(funs(sum(!is.na(.))))
head(A8summ_coef)

A8_SUM <- cbind(A8summ_coef[,1:2],
  exva = rowMeans(A8summ_coef[,4]),
  novar = rowMeans(A8summ_coef[,5]))

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
  )  + ggtitle("s=0.8")

A8BIC_PLOT
ggarrange(A05BIC_PLOT, A1BIC_PLOT, A3BIC_PLOT, A8BIC_PLOT, ncol = 2, nrow = 2) 

