g_deviance <- function(nu, data_in) {
  logLik(glm(y~.,family = binomial(Gosset(nu)),data = data_in, maxit = 100))
}

num_exp <- function(n_obs,s,n_out,df_stu,seeed){
  set.seed(seeed)
  
  n_obs <- n_obs*1.3
  alpha <- 0
  beta1 <- 1
  Sigma <- matrix(0, nrow = 1, ncol = 1)
  diag(Sigma) <- 1
  Ex <- cbind(MASS::mvrnorm(n_obs, rep(0, 1), Sigma))
  eta1 <- alpha + as.matrix(Ex) %*% (beta1)
  h_x_beta <- S(eta1,s)
  y <- rbinom(n = n_obs,size = 1,prob = h_x_beta)  
  dat2 <- data.frame(y,Ex)
  dat2$y <- as.factor(dat2$y)
  dat3 <- dat2
  
  n_out_f <- n_out
  Sigma2 <- matrix(0, nrow = 1, ncol = 1)
  diag(Sigma2) <- 0.5
  if(n_out_f > 0){
    if(n_out_f==1){
      outliers <- data.frame(y=rep(0,n_out_f),Ex=data.frame(MASS::mvrnorm(n = 2, rep(2, 1), Sigma2))) 
      colnames(outliers) <- colnames(dat2)
      dat3 <- rbind(dat2,outliers[1,])
    }else{
      outliers <- data.frame(y=rep(0,n_out_f),Ex=data.frame(MASS::mvrnorm(n = n_out_f, rep(2, 1), Sigma2))) 
      colnames(outliers) <- colnames(dat2)
      dat3 <- rbind(dat2,outliers)
    }
  }
  
  dat_test <- dat3[1:30,]
  n_test <- nrow(dat_test)
  dat3 <- dat3[31:nrow(dat3),]
  
  # nu_expl <- c(seq(0.2,8,0.05))
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
  
  mod_log <- glm(formula = y ~ . , data = dat3, family = binomial(link = "logit"), maxit = 100)
  
  log_nu3 <- tryCatch(
    glm(formula = y ~ ., data = dat3, family = binomial(link = Gosset(df_stu))),
    warning = 
      function(e){
        return(0)
      })
  # No models without convergence
  if(length(log_nu3) == 1 ){stop("")}
  
  mod_stu <- glm(formula = y ~ ., data = dat3, family = binomial(link = Gosset(df_stu)), maxit = 100)
  
  pred_log <- ifelse(predict(mod_log, dat_test, type = "response")<0.5,0,1)
  tab_log <- table(dat_test$y,pred_log)
  acc_log <- sum(diag(tab_log))/n_test
  
  pred_stu <- ifelse(predict(mod_stu, dat_test, type = "response")<0.5,0,1)
  tab_stu <- table(dat_test$y,pred_stu)
  acc_stu <- sum(diag(tab_stu))/n_test
  
  
  return(list(
    log_logit = logLik(mod_log),
    log_stu = logLik(mod_stu),
    grille = log_nu,
    df_stu = df_stu,
    acc_log = acc_log,
    acc_stu = acc_stu
  ))
}


# 0.05 --------------------------------------------------------------------

num_exp(n_obs = 100, s = 0.05, n_out = 2, df_stu = 0.8, seeed = 1)

nsim <- 100
A005log_logit <- A005log_stu <- A005grille <- A005df_stu <- A005acc_log <- A005acc_stu <- NULL
variable <- 100
for (outlier in 0:20) {
  # for (outlier in seq(0,20,2)) {
  # outlier <- 0
  for (nsim_q in 1:nsim) {
    # variable <- 1
    variable <- variable+1
    print(paste0(outlier,"A",variable))
    mtry <- try(num_exp(n_obs = 100, s = 0.05, n_out = outlier, df_stu = 0.8, seeed = variable))
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp(n_obs = 100, s = 0.05, n_out = outlier, df_stu = 0.8, seeed = variable))
    } 
    mode1 <- mtry
    # mode1 <- num_exp(n_obs = 100, s = 0.1, n_out = outlier, df_stu = 0.8, seeed = variable)
    A005log_logit <- rbind(A005log_logit,cbind(mode1[["log_logit"]],outlier,variable))
    A005log_stu <- rbind(A005log_stu,cbind(mode1[["log_stu"]],outlier,variable))
    A005grille <- rbind(A005grille,cbind(t(mode1[["grille"]]),outlier,variable))
    A005df_stu <- rbind(A005df_stu,cbind(mode1[["df_stu"]],outlier,variable))
    A005acc_log <- rbind(A005acc_log,cbind(mode1[["acc_log"]],outlier,variable))
    A005acc_stu <- rbind(A005acc_stu,cbind(mode1[["acc_stu"]],outlier,variable))
  }
  # print(outlier)
}

View(A005df_stu)

# Plot Loglikelihoods -----------------------------------------------------
A005log_logit <- as.data.frame(A005log_logit)
colnames(A005log_logit)[1] <- "loglik"
A005log_logit$link <- "logit"

A005log_stu <- as.data.frame(A005log_stu)
colnames(A005log_stu)[1] <- "loglik"
A005log_stu$link <- "stu"
# A005log_logit <- A005log_stu

A005logliks <- rbind(A005log_logit,A005log_stu)
str(A005logliks)

A005loglik1_plot <- ggplot(data = A005logliks[A005logliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("Log-likelihood") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(strip.background = element_rect(fill="yellow"),
    panel.background = element_rect(fill = "yellow", colour = "grey50"),
    plot.background = element_rect(fill = "yellow", colour = "yellow"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "yellow"),
    legend.key = element_rect(fill = NA),
    legend.position = c(0.9,0.8)
  )  + ggtitle("$d=0.05$")

A005loglik1_plot


# Plot Accuracies ---------------------------------------------------------

A005acc_logit <- as.data.frame(A005acc_log)
colnames(A005acc_logit)[1] <- "loglik"
A005acc_logit$link <- "logit"

A005acc_stu <- as.data.frame(A005acc_stu)
colnames(A005acc_stu)[1] <- "loglik"
A005acc_stu$link <- "stu"
# A005acc_logit <- A005acc_stu

A005accliks <- rbind(A005acc_logit,A005acc_stu)

A005acc_plot <- ggplot(data = A005accliks[A005accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("Accuracy") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.05$")

A005acc_plot
# 0.1 --------------------------------------------------------------------

num_exp(n_obs = 100, s = 0.1, n_out = 2, df_stu = 0.8, seeed = 1)

nsim <- 100
A01log_logit <- A01log_stu <- A01grille <- A01df_stu <- A01acc_log <- A01acc_stu <- NULL
variable <- 100
for (outlier in 0:20) {
  # for (outlier in seq(0,20,2)) {
  # outlier <- 0
  for (nsim_q in 1:nsim) {
    # variable <- 1
    variable <- variable+1
    print(paste0(outlier,"A",variable))
    mtry <- try(num_exp(n_obs = 100, s = 0.1, n_out = outlier, df_stu = 0.8, seeed = variable))
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp(n_obs = 100, s = 0.1, n_out = outlier, df_stu = 0.8, seeed = variable))
    } 
    mode1 <- mtry
    # mode1 <- num_exp(n_obs = 100, s = 0.1, n_out = outlier, df_stu = 0.8, seeed = variable)
    A01log_logit <- rbind(A01log_logit,cbind(mode1[["log_logit"]],outlier,variable))
    A01log_stu <- rbind(A01log_stu,cbind(mode1[["log_stu"]],outlier,variable))
    A01grille <- rbind(A01grille,cbind(t(mode1[["grille"]]),outlier,variable))
    A01df_stu <- rbind(A01df_stu,cbind(mode1[["df_stu"]],outlier,variable))
    A01acc_log <- rbind(A01acc_log,cbind(mode1[["acc_log"]],outlier,variable))
    A01acc_stu <- rbind(A01acc_stu,cbind(mode1[["acc_stu"]],outlier,variable))
  }
  # print(outlier)
}

View(A01df_stu)

# Plot Loglikelihoods -----------------------------------------------------
A01log_logit <- as.data.frame(A01log_logit)
colnames(A01log_logit)[1] <- "loglik"
A01log_logit$link <- "logit"

A01log_stu <- as.data.frame(A01log_stu)
colnames(A01log_stu)[1] <- "loglik"
A01log_stu$link <- "stu"
# A01log_logit <- A01log_stu

A01logliks <- rbind(A01log_logit,A01log_stu)
str(A01logliks)

A01loglik1_plot <- ggplot(data = A01logliks[A01logliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("Log-likelihood") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(strip.background = element_rect(fill="yellow"),
    panel.background = element_rect(fill = "yellow", colour = "grey50"),
    plot.background = element_rect(fill = "yellow", colour = "yellow"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "yellow"),
    legend.key = element_rect(fill = NA),
    legend.position = c(0.9,0.8)
  )  + ggtitle("$d=0.1$")

A01loglik1_plot


# Plot Accuracies ---------------------------------------------------------

A01acc_logit <- as.data.frame(A01acc_log)
colnames(A01acc_logit)[1] <- "loglik"
A01acc_logit$link <- "logit"

A01acc_stu <- as.data.frame(A01acc_stu)
colnames(A01acc_stu)[1] <- "loglik"
A01acc_stu$link <- "stu"
# A01acc_logit <- A01acc_stu

A01accliks <- rbind(A01acc_logit,A01acc_stu)

A01acc_plot <- ggplot(data = A01accliks[A01accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("Accuracy") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.1$")

A01acc_plot

# 0.3 --------------------------------------------------------------------

num_exp(n_obs = 100, s = 0.3, n_out = 2, df_stu = 0.8, seeed = 1)

nsim <- 100
A03log_logit <- A03log_stu <- A03grille <- A03df_stu <- A03acc_log <- A03acc_stu <- NULL
variable <- 100
for (outlier in 0:20) {
  # for (outlier in seq(0,20,2)) {
  # outlier <- 0
  for (nsim_q in 1:nsim) {
    # variable <- 1
    variable <- variable+1
    print(paste0(outlier,"A",variable))
    mtry <- try(num_exp(n_obs = 100, s = 0.3, n_out = outlier, df_stu = 0.8, seeed = variable))
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp(n_obs = 100, s = 0.3, n_out = outlier, df_stu = 0.8, seeed = variable))
    } 
    mode1 <- mtry
    # mode1 <- num_exp(n_obs = 100, s = 0.3, n_out = outlier, df_stu = 0.8, seeed = variable)
    A03log_logit <- rbind(A03log_logit,cbind(mode1[["log_logit"]],outlier,variable))
    A03log_stu <- rbind(A03log_stu,cbind(mode1[["log_stu"]],outlier,variable))
    A03grille <- rbind(A03grille,cbind(t(mode1[["grille"]]),outlier,variable))
    A03df_stu <- rbind(A03df_stu,cbind(mode1[["df_stu"]],outlier,variable))
    A03acc_log <- rbind(A03acc_log,cbind(mode1[["acc_log"]],outlier,variable))
    A03acc_stu <- rbind(A03acc_stu,cbind(mode1[["acc_stu"]],outlier,variable))
  }
  # print(outlier)
}

View(A03df_stu)

# Plot Loglikelihoods -----------------------------------------------------
A03log_logit <- as.data.frame(A03log_logit)
colnames(A03log_logit)[1] <- "loglik"
A03log_logit$link <- "logit"

A03log_stu <- as.data.frame(A03log_stu)
colnames(A03log_stu)[1] <- "loglik"
A03log_stu$link <- "stu"
# A03log_logit <- A03log_stu

A03logliks <- rbind(A03log_logit,A03log_stu)
str(A03logliks)

A03loglik1_plot <- ggplot(data = A03logliks[A03logliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("Log-likelihood") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(strip.background = element_rect(fill="yellow"),
    panel.background = element_rect(fill = "yellow", colour = "grey50"),
    plot.background = element_rect(fill = "yellow", colour = "yellow"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "yellow"),
    legend.key = element_rect(fill = NA),
    legend.position = c(0.9,0.8)
  )  + ggtitle("$d=0.3$")

A03loglik1_plot


# Plot Accuracies ---------------------------------------------------------

A03acc_logit <- as.data.frame(A03acc_log)
colnames(A03acc_logit)[1] <- "loglik"
A03acc_logit$link <- "logit"

A03acc_stu <- as.data.frame(A03acc_stu)
colnames(A03acc_stu)[1] <- "loglik"
A03acc_stu$link <- "stu"
# A03acc_logit <- A03acc_stu

A03accliks <- rbind(A03acc_logit,A03acc_stu)

A03acc_plot <- ggplot(data = A03accliks[A03accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("Accuracy") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(strip.background = element_rect(fill="yellow"),
    panel.background = element_rect(fill = "yellow", colour = "grey50"),
    plot.background = element_rect(fill = "yellow", colour = "yellow"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )  + ggtitle("$d=0.3$")

A03acc_plot

# 0.8 --------------------------------------------------------------------

num_exp(n_obs = 100, s = 0.8, n_out = 2, df_stu = 0.8, seeed = 1)

nsim <- 100
A08log_logit <- A08log_stu <- A08grille <- A08df_stu <- A08acc_log <- A08acc_stu <- NULL
variable <- 100
for (outlier in 0:20) {
  # for (outlier in seq(0,20,2)) {
  # outlier <- 0
  for (nsim_q in 1:nsim) {
    # variable <- 1
    variable <- variable+1
    print(paste0(outlier,"A",variable))
    mtry <- try(num_exp(n_obs = 100, s = 0.8, n_out = outlier, df_stu = 0.8, seeed = variable))
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      mtry <- try(num_exp(n_obs = 100, s = 0.8, n_out = outlier, df_stu = 0.8, seeed = variable))
    } 
    mode1 <- mtry
    # mode1 <- num_exp(n_obs = 100, s = 0.8, n_out = outlier, df_stu = 0.8, seeed = variable)
    A08log_logit <- rbind(A08log_logit,cbind(mode1[["log_logit"]],outlier,variable))
    A08log_stu <- rbind(A08log_stu,cbind(mode1[["log_stu"]],outlier,variable))
    A08grille <- rbind(A08grille,cbind(t(mode1[["grille"]]),outlier,variable))
    A08df_stu <- rbind(A08df_stu,cbind(mode1[["df_stu"]],outlier,variable))
    A08acc_log <- rbind(A08acc_log,cbind(mode1[["acc_log"]],outlier,variable))
    A08acc_stu <- rbind(A08acc_stu,cbind(mode1[["acc_stu"]],outlier,variable))
  }
  # print(outlier)
}

View(A08df_stu)


# Plot Loglikelihoods -----------------------------------------------------
A08log_logit <- as.data.frame(A08log_logit)
colnames(A08log_logit)[1] <- "loglik"
A08log_logit$link <- "logit"

A08log_stu <- as.data.frame(A08log_stu)
colnames(A08log_stu)[1] <- "loglik"
A08log_stu$link <- "stu"
# A08log_logit <- A08log_stu

A08logliks <- rbind(A08log_logit,A08log_stu)
str(A08logliks)

A08loglik1_plot <- ggplot(data = A08logliks[A08logliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("Log-likelihood") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(strip.background = element_rect(fill="yellow"),
    panel.background = element_rect(fill = "yellow", colour = "grey50"),
    plot.background = element_rect(fill = "yellow", colour = "yellow"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "yellow"),
    legend.key = element_rect(fill = NA),
    legend.position = c(0.9,0.8)
  )  + ggtitle("$d=0.8$")

A08loglik1_plot


# Plot Accuracies ---------------------------------------------------------

A08acc_logit <- as.data.frame(A08acc_log)
colnames(A08acc_logit)[1] <- "loglik"
A08acc_logit$link <- "logit"

A08acc_stu <- as.data.frame(A08acc_stu)
colnames(A08acc_stu)[1] <- "loglik"
A08acc_stu$link <- "stu"
# A08acc_logit <- A08acc_stu

A08accliks <- rbind(A08acc_logit,A08acc_stu)

A08acc_plot <- ggplot(data = A08accliks[A08accliks$outlier %in% seq(0,20,2),], aes(x=factor(outlier),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("Accuracy") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  )  + ggtitle("$d=0.8$")

A08acc_plot


ggarrange(A005loglik1_plot, A01loglik1_plot, A03loglik1_plot, A08loglik1_plot, ncol = 2, nrow = 2) 
ggarrange(A005acc_plot, A01acc_plot, A03acc_plot, A08acc_plot, ncol = 2, nrow = 2) 

tikz(file = "tex/Simu_outliers_log_lik.tex", standAlone=F, width = 7.5, height = 4.5)
endoffile <- dev.off()

tikz(file = "tex/Simu_outliers_acc.tex", standAlone=F, width = 7.5, height = 4.5)
endoffile <- dev.off()

tikz(file = "tex/Simu_outliers_log_lik_005.tex", standAlone=F, width = 5.5, height = 2.8)
endoffile <- dev.off()

tikz(file = "tex/Simu_outliers_log_lik_001.tex", standAlone=F, width = 5.5, height = 2.8)
endoffile <- dev.off()

tikz(file = "tex/Simu_outliers_log_lik_03.tex", standAlone=F, width = 5.5, height = 2.8)
endoffile <- dev.off()

tikz(file = "tex/Simu_outliers_log_lik_08.tex", standAlone=F, width = 5.5, height = 2.8)
endoffile <- dev.off()

# Plot outliers -----------------------------------------------------------

sum(A005acc_log<A005acc_stu)/210

View(A005grille)
colnames(A005grille) <- c(nu_expl,"outlier", "seed")
dat_plot <- cbind(A005grille,A005df_stu)
dat_plot_8 <- dat_plot
# dat_plot_8 <- dat_plot_8[dat_plot_8[,ncol(dat_plot)-1] == 8,]
View(dat_plot_8)
colnames(dat_plot_8)
head(dat_plot_8)
dat_plot_8 <- as.data.frame(dat_plot_8)
dat_plot_8 <- dat_plot_8[,-ncol(dat_plot_8)]
dat_plot_8 <- dat_plot_8[,-ncol(dat_plot_8)]
dat_plot_8 <- dat_plot_8[,-ncol(dat_plot_8)]
dat_plot_8 <- dat_plot_8[,-ncol(dat_plot_8)]
# dat_plot_8 <- dat_plot_8[,-ncol(dat_plot_8)]
dat_plot_8$Ind <- 1:nrow(dat_plot_8)

library(reshape2)

long_8 <- melt(dat_plot_8, id.vars=c("outlier", "Ind"))

str(long_8)
long_8$variable <- as.numeric(as.character(long_8$variable))
head(long_8)

set.seed(9)
dat_r <- long_8 %>% 
  group_by(outlier) %>%
  do(sample_n(.,2))
dat_r$Ind

ggplot(data = long_8,mapping = aes(x = variable, y = value, group = factor(Ind), color = outlier))+
  # ggplot(data = long_8[long_8$Ind %in% dat_r$Ind,],mapping = aes(x = variable, y = value, group = factor(Ind), color = outlier))+
  # geom_point()+
  geom_line() +
  xlim(0.2,8)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  xlab("v") + ylab("Log-likelihood") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white")
    # legend.position = "none"
  )  
# + ggtitle("s=0.05")

ggplot(data = long_8[long_8$outlier<15,],mapping = aes(x = variable, y = value, group = factor(Ind), color = outlier))+
  # geom_point()+
  geom_line() +
  xlab("nu") + ylab("Loglikelihood") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white")
    # legend.position = "none"
  )  + 
  xlim(0.3,4)










