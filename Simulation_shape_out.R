
# Functions optimization --------------------------------------------------

g_deviance <- function(nu, data_in) {
  logLik(glm(response~.,family = binomial(Gosset(nu)),data = data_in, maxit = 100))
}


# Simulation cases --------------------------------------------------------

shape1  <- 2
shape2 <- 2

beta1 <- 1
alpha <- 0
amplitude_.5 <- 1
n_obs <- 250
n_test <- 250*0.2
n_training <- 250*0.8
p_var <- 1

z_beta <- rbeta(n = n_obs, shape1, shape2)

ggplot(data.frame(z_beta), aes(x = z_beta )) +
  # geom_histogram()+
  geom_density()+
  # geom_point(y = 0)+
  geom_jitter(aes(y=0),width = 0.0000000000005, height = 8)

Sigma <- matrix(0.5, nrow = p_var, ncol = p_var)
diag(Sigma) <- 1
Ex1 <- cbind(mvrnorm(n_obs, rep(0, p_var), Sigma))

Ex <- ((2*amplitude_.5*z_beta) / (beta1*(max(z_beta)-min(z_beta)))) + ((-amplitude_.5-alpha) / beta1)
Ex <- cbind(Ex,Ex1[,-1])

beta1 <- c(1, rep(0, p_var - 1))

eta1 <- alpha + as.matrix(Ex) %*% (beta1)

S<-function(x,s){ return  (1/(1+exp(-x/s)))}
s <-  0.05

h_x_beta <- S(eta1,s)

y1 <- rbinom(n = n_obs,size = 1,prob = h_x_beta)    

data_sim <- as.data.frame(cbind(Ex, response = y1))

ggplot(data = data_sim) +
  # geom_point(aes(x = z_beta, y = y1, color = factor(y1)))
  geom_jitter(aes(x = Ex, y = response, color = factor(response)),width = 0.0000000000005, height = 0.05)

# One Simulation  ---------------------------------------------------------


shape1  <- 2
shape2 <- 2

beta1 <- 1
alpha <- 0
amplitude_.5 <- 5
n_obs <- 250
n_test <- 250*0.2
n_training <- 250*0.8
p_var <- 1



s <-  0.05

n_sim <- 10

df_stu <- loglik_logit <- loglik_stu <- accu_logit <- accu_stu <- matrix(nrow = n_sim, ncol = 1)


simu_f <- function(p_out,s,n_sim){
  # p_out <- 5
  # s <- 0.5
  n_out <- n_training*p_out/100
  sim <- 0
  
  while (sim < n_sim) {
    # print(sim)
    sim <- sim+1
    z_beta <- rbeta(n = n_obs, shape1, shape2)
    Sigma <- matrix(0.5, nrow = p_var, ncol = p_var)
    diag(Sigma) <- 1
    # Ex1 <- cbind(mvrnorm(n_obs, rep(0, p_var), Sigma))
    Ex <- (amplitude_.5-alpha)*((2*z_beta) - 1) / beta1
    # Ex <- ((2*amplitude_.5*z_beta) / (beta1*(max(z_beta)-min(z_beta)))) + ((-amplitude_.5-alpha) / beta1)
    # Ex <- cbind(Ex,Ex1[,-1])
    
    beta1 <- c(1, rep(0, p_var - 1))
    eta1 <- alpha + as.matrix(Ex) %*% (beta1)
    
    
    h_x_beta <- S(eta1,s)
    y1 <- rbinom(n = n_obs,size = 1,prob = h_x_beta)    
    
    data_sim <- as.data.frame(cbind(Ex, response = y1))
    dat_training <- data_sim[1:n_training,]
    dat_test <- data_sim[(n_training+1):n_obs,]
    
    if (n_out > 0) {
      dat_training[
        rownames(dat_training[dat_training$Ex %in% sort(dat_training[,"Ex"], F)[1:n_out],])
        ,"response"] <- rep(1,n_out)
    }
    
    result_nu <- optimize(g_deviance, c(0.25, 4), maximum = T, dat_training)
    df_stu[sim] <- result_nu$maximum
    
    result_nu <- tryCatch(glm(response~.,family = binomial(Gosset(result_nu$maximum)),data = dat_training), 
                          warning = 
                            function(e){
                              return(0)
                            })
    
    if(length(result_nu) == 1 ){sim <- sim-1}
    # print(sim)
    if(length(result_nu) == 1 ){next}
    
    # mod_stu <- glm(response~.,family = binomial(Gosset(result_nu$maximum)),data = dat_training)
    mod_stu <- result_nu
    loglik_stu[sim] <- logLik(mod_stu)
    pred_stu <- ifelse(predict(mod_stu, dat_test, type = "response")<0.5,0,1)
    tab_stu <- table(dat_test$response,pred_stu)
    acc_stu <- sum(diag(tab_stu))/n_test
    accu_stu[sim] <- acc_stu
    
    mod_logit <- glm(response~.,family = binomial("logit"),data = dat_training, maxit = 100)
    loglik_logit[sim] <- logLik(mod_logit)
    pred_logit <- ifelse(predict(mod_logit, dat_test, type = "response")<0.5,0,1)
    tab_logit <- table(dat_test$response,pred_logit)
    acc_logit <- sum(diag(tab_logit))/n_test
    accu_logit[sim] <- acc_logit  
    
  }
  
  log_lik_res <- rbind(data.frame(loglik=loglik_logit,link="logistic",pi=s,n_out =p_out),
                       data.frame(loglik=loglik_stu,link="student",pi=s,n_out =p_out))
  
  accu_res <- rbind(data.frame(accuracy=accu_logit,link="logistic",pi=s,n_out =p_out),
                    data.frame(accuracy=accu_stu,link="student",pi=s,n_out =p_out))
  
  df_res <- rbind(data.frame(df=df_stu,pi=s,n_out =p_out))
  
  return(list(log_lik_res,accu_res,df_res))
}

dat_training$s <- 0.5
dat_training1 <- dat_training
dat_training2 <- dat_training
dat_training3 <- dat_training

dat_training_c <- rbind(dat_training1, dat_training2, dat_training3)

ggplot(data = dat_training_c) +
  # geom_point(aes(x = z_beta, y = y1, color = factor(y1)))
  geom_jitter(aes(x = Ex, y = response, color = factor(response)),width = 0.0000000000005, height = 0.05)+
  xlab("X1") + ylab("y") +
  scale_y_continuous(breaks = c(0,1),
                     labels = c("0","1"))+
  scale_color_manual("", labels = c("0", "1"), values = c("#4E84C4", "#52854C"))+
  facet_grid(cols = vars(s))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none"
  ) 

# tikz(file = "Scatter_plot_ex.tex", standAlone=F,width = 7, height = 3)
# endoffile <- dev.off() 


# s1 ----------------------------------------------------------------------

log_lik_p_s_t <- data.frame(matrix(NA, nrow = 1, ncol = 4))
colnames(log_lik_p_s_t) <- c("loglik", "link", "pi", "n_out")
acc_p_s_t <- matrix(NA, nrow = 1, ncol = 4)
colnames(acc_p_s_t) <- c("accuracy", "link", "pi", "n_out")
df_lik_p_s_t <- matrix(NA, nrow = 1, ncol = 3)
colnames(df_lik_p_s_t) <- c("df", "pi", "n_out")


for (n_outliers in seq(0,20,2)) {
  # n_outliers <- 0
  p_s <- simu_f(p_out = n_outliers, 0.05, 100)
  
  log_lik_p_s <- p_s[[1]]
  acc_p_s <- p_s[[2]]
  df_lik_p_s <- p_s[[3]]
  
  log_lik_p_s_t <- rbind(log_lik_p_s_t,log_lik_p_s)
  acc_p_s_t <- rbind(acc_p_s_t,acc_p_s)
  df_lik_p_s_t <- rbind(df_lik_p_s_t,df_lik_p_s)
  
  print(n_outliers)
}

log_lik_p_s_t <- log_lik_p_s_t[-1,]
acc_p_s_t <- acc_p_s_t[-1,]
df_lik_p_s_t <- df_lik_p_s_t[-1,]

head(log_lik_p_s_t)
summary(log_lik_p_s_t)

loglik1_plot <- ggplot(data = log_lik_p_s_t, aes(x=factor(n_out),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("Log-Likelihood") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = c(0.88,0.8)
  )  + ggtitle("Log-likelihood profiles")

loglik1_plot

tikz(file = "Log_Lik_s11.tex", standAlone=F,width = 6.5, height = 3)
endoffile <- dev.off() 

acc1_plot <- ggplot(data = acc_p_s_t, aes(x=factor(n_out),y=accuracy,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("Accuracy") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        # legend.position = c(0.88,0.8),
        legend.position = "none"
  ) + ggtitle("Accuracy")



head(df_lik_p_s_t)

out1_plot <- ggplot(data = df_lik_p_s_t, aes(x=factor(n_out),y=df))+
  geom_jitter(color = "#00AFBB", height = 0.1)+
  # facet_grid(cols = vars(n_out))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none"
  ) + ggtitle("Selected nu") + xlab("of outliers")

out1_plot

library(ggpubr)



ggarrange(loglik1_plot,                                                 # First row with scatter plot
          ggarrange(acc1_plot, out1_plot, ncol = 2, nrow = 1
                    # labels = c("Accuracy", "Selected nu")
          ), # Second row with box and dot plots
          nrow = 2
          # labels = "Log-likelihood profiles"                                        # Labels of the scatter plot
) 

tikz(file = "S1_ll_acc_df_ALL.tex", standAlone=F,width = 7, height = 8)
endoffile <- dev.off() 

# s2 ----------------------------------------------------------------------

log_lik_p_s_t_2 <- data.frame(matrix(NA, nrow = 1, ncol = 4))
colnames(log_lik_p_s_t_2) <- c("loglik", "link", "pi", "n_out")
acc_p_s_t_2 <- matrix(NA, nrow = 1, ncol = 4)
colnames(acc_p_s_t_2) <- c("accuracy", "link", "pi", "n_out")
df_lik_p_s_t_2 <- matrix(NA, nrow = 1, ncol = 3)
colnames(df_lik_p_s_t_2) <- c("df", "pi", "n_out")

for (n_outliers in seq(0,20,2)) {
  # n_outliers <- 0
  p_s <- simu_f(p_out = n_outliers, 0.1, 100)
  
  log_lik_p_s <- p_s[[1]]
  acc_p_s <- p_s[[2]]
  df_lik_p_s <- p_s[[3]]
  
  log_lik_p_s_t_2 <- rbind(log_lik_p_s_t_2,log_lik_p_s)
  acc_p_s_t_2 <- rbind(acc_p_s_t_2,acc_p_s)
  df_lik_p_s_t_2 <- rbind(df_lik_p_s_t_2,df_lik_p_s)
  
  print(n_outliers)
}

log_lik_p_s_t_2 <- log_lik_p_s_t_2[-1,]
acc_p_s_t_2 <- acc_p_s_t_2[-1,]
df_lik_p_s_t_2 <- df_lik_p_s_t_2[-1,]

head(log_lik_p_s_t)
summary(log_lik_p_s_t)

ggplot(data = log_lik_p_s_t_2, aes(x=factor(n_out),y=loglik,fill=link)) +
  geom_boxplot() +
  xlab("of outliers") + ylab("Log-Likelihood") +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = c(0.88,0.8)
  ) 

tikz(file = "Log_Lik_s2.tex", standAlone=F,width = 6.5, height = 3)
endoffile <- dev.off() 

ggplot(data = acc_p_s_t_2, aes(x=factor(n_out),y=accuracy,fill=link))+
  geom_boxplot()+
  # facet_grid(cols = vars(n_out))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white")
  ) 

# s3 ----------------------------------------------------------------------

log_lik_p_s_t_3 <- data.frame(matrix(NA, nrow = 1, ncol = 4))
colnames(log_lik_p_s_t_3) <- c("loglik", "link", "pi", "n_out")
acc_p_s_t_3 <- matrix(NA, nrow = 1, ncol = 4)
colnames(acc_p_s_t_3) <- c("accuracy", "link", "pi", "n_out")
df_lik_p_s_t_3 <- matrix(NA, nrow = 1, ncol = 3)
colnames(df_lik_p_s_t_3) <- c("df", "pi", "n_out")

for (n_outliers in seq(0,20,2)) {
  # n_outliers <- 0
  p_s <- simu_f(p_out = n_outliers, 0.3, 100)
  
  log_lik_p_s <- p_s[[1]]
  acc_p_s <- p_s[[2]]
  df_lik_p_s <- p_s[[3]]
  
  log_lik_p_s_t_3 <- rbind(log_lik_p_s_t_3,log_lik_p_s)
  acc_p_s_t_3 <- rbind(acc_p_s_t_3,acc_p_s)
  df_lik_p_s_t_3 <- rbind(df_lik_p_s_t_3,df_lik_p_s)
  
  print(n_outliers)
}

log_lik_p_s_t_3 <- log_lik_p_s_t_3[-1,]
acc_p_s_t_3 <- acc_p_s_t_3[-1,]
df_lik_p_s_t_3 <- df_lik_p_s_t_3[-1,]

head(log_lik_p_s_t)
summary(log_lik_p_s_t)

# s4 ----------------------------------------------------------------------

log_lik_p_s_t_4 <- data.frame(matrix(NA, nrow = 1, ncol = 4))
colnames(log_lik_p_s_t_4) <- c("loglik", "link", "pi", "n_out")
acc_p_s_t_4 <- matrix(NA, nrow = 1, ncol = 4)
colnames(acc_p_s_t_4) <- c("accuracy", "link", "pi", "n_out")
df_lik_p_s_t_4 <- matrix(NA, nrow = 1, ncol = 3)
colnames(df_lik_p_s_t_4) <- c("df", "pi", "n_out")

for (n_outliers in seq(0,20,2)) {
  # n_outliers <- 0
  p_s <- simu_f(p_out = n_outliers, 0.8, 100)
  
  log_lik_p_s <- p_s[[1]]
  acc_p_s <- p_s[[2]]
  df_lik_p_s <- p_s[[3]]
  
  log_lik_p_s_t_4 <- rbind(log_lik_p_s_t_4,log_lik_p_s)
  acc_p_s_t_4 <- rbind(acc_p_s_t_4,acc_p_s)
  df_lik_p_s_t_4 <- rbind(df_lik_p_s_t_4,df_lik_p_s)
  
  print(n_outliers)
}

log_lik_p_s_t_4 <- log_lik_p_s_t_4[-1,]
acc_p_s_t_4 <- acc_p_s_t_4[-1,]
df_lik_p_s_t_4 <- df_lik_p_s_t_4[-1,]


# graficos ----------------------------------------------------------------

loglik1_plot
loglik2_plot
loglik3_plot
loglik4_plot

ggarrange(loglik1_plot,
          loglik2_plot,
          loglik3_plot,
          loglik4_plot,
          nrow = 2, ncol = 2) 

tikz(file = "Log_Lik.tex", standAlone=F,width = 6.5, height = 4.2)
endoffile <- dev.off() 

loglik1_plot <- ggplot(data = log_lik_p_s_t, aes(x=factor(n_out),y=loglik,fill=link)) +
  geom_boxplot() +
  # xlab("") + ylab("") +
  scale_y_continuous( name = NULL)+
  scale_x_discrete(name = NULL)+
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none"
  ) + ggtitle("s=0.05")

loglik2_plot <- ggplot(data = log_lik_p_s_t_2, aes(x=factor(n_out),y=loglik,fill=link)) +
  geom_boxplot() +
  # xlab("") + ylab("") +
  scale_y_continuous( name = NULL)+
  scale_x_discrete(name = NULL)+
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none"
  ) + ggtitle("s=0.1")
loglik3_plot <- ggplot(data = log_lik_p_s_t_3, aes(x=factor(n_out),y=loglik,fill=link)) +
  geom_boxplot() +
  # xlab("") + ylab("") +
  scale_y_continuous( name = NULL)+
  scale_x_discrete(name = NULL)+
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none"
  ) + ggtitle("s=0.3")
loglik4_plot <- ggplot(data = log_lik_p_s_t_4, aes(x=factor(n_out),y=loglik,fill=link)) +
  geom_boxplot() +
  # xlab("") + ylab("") +
  scale_y_continuous( name = NULL)+
  scale_x_discrete(name = NULL)+
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none"
  ) + ggtitle("s=0.8")
# + ggtitle("Log-likelihood profiles")

tikz(file = "Log_Lik_s11.tex", standAlone=F,width = 6.5, height = 3)
endoffile <- dev.off() 


# Nu ----------------------------------------------------------------------

out1_plot <- ggplot(data = df_lik_p_s_t, aes(x=factor(n_out),y=df))+
  geom_jitter(color = "#00AFBB", height = 0.1)+
  # facet_grid(cols = vars(n_out))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none"
  ) +
  scale_y_continuous( name = NULL)+
  scale_x_discrete(name = NULL)+
  ggtitle("s=0.05")

out1_plot

out2_plot <- ggplot(data = df_lik_p_s_t_2, aes(x=factor(n_out),y=df))+
  geom_jitter(color = "#00AFBB", height = 0.1)+
  # facet_grid(cols = vars(n_out))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none"
  ) +
  scale_y_continuous( name = NULL)+
  scale_x_discrete(name = NULL)+
  ggtitle("s=0.1")

out2_plot

out3_plot <- ggplot(data = df_lik_p_s_t_3, aes(x=factor(n_out),y=df))+
  geom_jitter(color = "#00AFBB", height = 0.1)+
  # facet_grid(cols = vars(n_out))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none"
  ) +
  scale_y_continuous( name = NULL)+
  scale_x_discrete(name = NULL)+
  ggtitle("s=0.3")

out4_plot <- ggplot(data = df_lik_p_s_t_4, aes(x=factor(n_out),y=df))+
  geom_jitter(color = "#00AFBB", height = 0.1)+
  # facet_grid(cols = vars(n_out))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none"
  ) +
  scale_y_continuous( name = NULL)+
  scale_x_discrete(name = NULL)+
  ggtitle("s=0.8")

ggarrange(out1_plot,
          out2_plot,
          out3_plot,
          out4_plot,
          nrow = 2, ncol = 2) 

tikz(file = "nu_selected.tex", standAlone=F,width = 6.5, height = 4.2)
endoffile <- dev.off() 

# Accuracies --------------------------------------------------------------

acc1_plot <- ggplot(data = acc_p_s_t, aes(x=factor(n_out),y=accuracy,fill=link)) +
  geom_boxplot() +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        # legend.position = c(0.88,0.8),
        legend.position = "none"
  ) +
  scale_y_continuous( name = NULL)+
  scale_x_discrete(name = NULL)+
  ggtitle("s=0.05")

acc2_plot <- ggplot(data = acc_p_s_t_2, aes(x=factor(n_out),y=accuracy,fill=link)) +
  geom_boxplot() +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        # legend.position = c(0.88,0.8),
        legend.position = "none"
  ) +
  scale_y_continuous( name = NULL)+
  scale_x_discrete(name = NULL)+
  ggtitle("s=0.1")

acc3_plot <- ggplot(data = acc_p_s_t_3, aes(x=factor(n_out),y=accuracy,fill=link)) +
  geom_boxplot() +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        # legend.position = c(0.88,0.8),
        legend.position = "none"
  ) +
  scale_y_continuous( name = NULL)+
  scale_x_discrete(name = NULL)+
  ggtitle("s=0.3")

acc4_plot <- ggplot(data = acc_p_s_t_4, aes(x=factor(n_out),y=accuracy,fill=link)) +
  geom_boxplot() +
  scale_fill_manual("", labels = c("$Logistic$", "$Student$"), values = c("#E69F00", "#00AFBB"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        # legend.position = c(0.88,0.8),
        legend.position = "none"
  ) +
  scale_y_continuous( name = NULL)+
  scale_x_discrete(name = NULL)+
  ggtitle("s=0.8")

ggarrange(acc1_plot,
          acc2_plot,
          acc3_plot,
          acc4_plot,
          nrow = 2, ncol = 2) 

tikz(file = "accuracies.tex", standAlone=F,width = 6.5, height = 4.2)
endoffile <- dev.off() 
