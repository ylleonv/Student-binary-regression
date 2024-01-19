nugridd <- seq(0.15,8,0.05)
num_exp <- function(n_obs,s,n_out,df_stu,seeed,out_y="meth2"){
  set.seed(seeed)
  alpha <- 0
  # beta1 <- c(1, rep(0, 1))
  beta1 <- 1
  Sigma <- matrix(0, nrow = 1, ncol = 1)
  diag(Sigma) <- 1
  Ex <- cbind(MASS::mvrnorm(n_obs, rep(0, 1), Sigma))
  eta1 <- alpha + as.matrix(Ex) %*% (beta1)
  h_x_beta <- pt(eta1,s)
  y <- rbinom(n = n_obs,size = 1,prob = h_x_beta)
  dat2 <- data.frame(y,Ex)
  dat2$y <- as.factor(dat2$y)
  dat3 <- dat2
  
  nu_expl <- nugridd
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
  
  # log_1 <- log_nu[nu_expl==1]
  # log_8 <- log_nu[nu_expl==8]
  
  # df_stu <- 8
  
  # if(log_1 > log_8){
  #   result_nu <- optimize(g_deviance, c(0.25, 1), maximum = T, dat3)
  #   df_stu_opt <- result_nu$maximum
  #   df_stu <- df_stu_opt
  # }
  
  return(list(
    log_nu = log_nu
    # df_stu = df_stu
  ))
}

A005coefs <- NULL
nsim <- 10
variable <- 20000
# for (s_stu in seq(0.5,4,0.3)) {
for (s_stu in c(0.3)) {
  # for (s_stu in seq(0.07,2,0.1)) {
  # outlier <- 0
  # s_stu <- 0.25
  # s_stu
  variable <- variable+1
  for (otro in 1:nsim) {
    # variable <- 1
    variable <- variable+1
    
    print(paste0(s_stu,"A",variable))
    
    mtry <- try(num_exp(n_obs = 100,s = s_stu,n_out = 0,df_stu = 0.8,seeed = variable))
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      print(paste0(s_stu,"A",variable))
      mtry <- try(num_exp(n_obs = 100,s = s_stu,n_out = 0,df_stu = 0.8,seeed = variable))
    } 
    
    mode1 <- mtry
    a2 <- data.frame(log_nu = mode1$log_nu)
    a2$grid <- nugridd
    a2$s_stu <- s_stu
    a2$variable <- variable
    A005coefs <- rbind(A005coefs,a2)
  }
}
A005coefs$v <- A005coefs$s_stu
plot03 <- ggplot(data = A005coefs,mapping = aes(x = grid, y = log_nu, group = factor(variable),
  color = factor(variable)
))+
  # ggplot(data = long_8[long_8$Ind %in% dat_r$Ind,],mapping = aes(x = variable, y = value, group = factor(Ind), color = outlier))+
  # geom_point()+
  geom_line() +
  # xlim(0.2,8)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  xlab("v") + ylab("Log-likelihood") + ggtitle("nu=0.3")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white")
    # legend.position = "none"
  ) 
plot03

# 0.3 ---------------------------------------------------------------------

A005coefs03 <- NULL
nsim <- 5
variable_i <- 0
variable_v <- c(10488,11856,10537,10869,20258)
# for (s_stu in seq(0.5,4,0.3)) {
for (s_stu in c(0.3)) {
  # variable <- variable+1
  for (otro in 1:nsim) {
    # variable <- 1
    variable_i <- variable_i+1
    variable <- variable_v[variable_i]
    
    print(paste0(s_stu,"A",variable))
    
    mtry <- try(num_exp(n_obs = 100,s = s_stu,n_out = 0,df_stu = 0.8,seeed = variable))
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      print(paste0(s_stu,"A",variable))
      mtry <- try(num_exp(n_obs = 100,s = s_stu,n_out = 0,df_stu = 0.8,seeed = variable))
    } 
    mode1 <- mtry
    a2 <- data.frame(log_nu = mode1$log_nu)
    a2$grid <- nugridd
    a2$s_stu <- s_stu
    a2$variable <- variable
    A005coefs03 <- rbind(A005coefs03,a2)
  }
}
A005coefs03$v <- A005coefs03$s_stu
head(A005coefs03)
plot03 <- ggplot(data = A005coefs03,mapping = aes(x = grid, y = log_nu, group = factor(variable)
  # color = factor(variable)
))+
  # ggplot(data = long_8[long_8$Ind %in% dat_r$Ind,],mapping = aes(x = variable, y = value, group = factor(Ind), color = outlier))+
  # geom_point()+
  geom_line() +
  # xlim(0.2,8)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  xlab("v") + 
  ylab("Log-likelihood") +
  ggtitle("v*=0.3")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white")
    # legend.position = "none"
  ) 
plot03

# 1.5 ---------------------------------------------------------------------

A005coefs15 <- NULL
nsim <- 6
variable_i <- 0
variable_v <- c(386,1558,1431,1356,10144,10080)
# for (s_stu in seq(0.5,4,0.3)) {
for (s_stu in c(1.5)) {
  # variable <- variable+1
  for (otro in 1:nsim) {
    # variable <- 1
    variable_i <- variable_i+1
    variable <- variable_v[variable_i]
    
    print(paste0(s_stu,"A",variable))
    
    mtry <- try(num_exp(n_obs = 100,s = s_stu,n_out = 0,df_stu = 0.8,seeed = variable))
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      print(paste0(s_stu,"A",variable))
      mtry <- try(num_exp(n_obs = 100,s = s_stu,n_out = 0,df_stu = 0.8,seeed = variable))
    } 
    mode1 <- mtry
    a2 <- data.frame(log_nu = mode1$log_nu)
    a2$grid <- nugridd
    a2$s_stu <- s_stu
    a2$variable <- variable
    A005coefs15 <- rbind(A005coefs15,a2)
  }
}
A005coefs15$v <- A005coefs15$s_stu
head(A005coefs15)
plot15 <- ggplot(data = A005coefs15,mapping = aes(x = grid, y = log_nu, group = factor(variable)
  # color = factor(variable)
))+
  # ggplot(data = long_8[long_8$Ind %in% dat_r$Ind,],mapping = aes(x = variable, y = value, group = factor(Ind), color = outlier))+
  # geom_point()+
  geom_line() +
  # xlim(0.2,8)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  xlab("v") + 
  ylab("Log-likelihood") +
  ggtitle("v*=1.5")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white")
    # legend.position = "none"
  ) 
plot15

# nu=6 --------------------------------------------------------------------

A005coefs6 <- NULL
nsim <- 5
variable_i <- 0
variable_v <- c(16,221,90,251,1013)
# for (s_stu in seq(0.5,4,0.3)) {
for (s_stu in c(6)) {
  # variable <- variable+1
  for (otro in 1:nsim) {
    # variable <- 1
    variable_i <- variable_i+1
    variable <- variable_v[variable_i]
    
    print(paste0(s_stu,"A",variable))
    
    mtry <- try(num_exp(n_obs = 100,s = s_stu,n_out = 0,df_stu = 0.8,seeed = variable))
    
    while (inherits(mtry, "try-error")) {
      variable <- variable+1
      print(paste0(s_stu,"A",variable))
      mtry <- try(num_exp(n_obs = 100,s = s_stu,n_out = 0,df_stu = 0.8,seeed = variable))
    } 
    mode1 <- mtry
    a2 <- data.frame(log_nu = mode1$log_nu)
    a2$grid <- nugridd
    a2$s_stu <- s_stu
    a2$variable <- variable
    A005coefs6 <- rbind(A005coefs6,a2)
  }
}
A005coefs6$v <- A005coefs6$s_stu
head(A005coefs6)
plot6 <- ggplot(data = A005coefs6,mapping = aes(x = grid, y = log_nu, group = factor(variable)
  # color = factor(variable)
))+
  # ggplot(data = long_8[long_8$Ind %in% dat_r$Ind,],mapping = aes(x = variable, y = value, group = factor(Ind), color = outlier))+
  # geom_point()+
  geom_line() +
  # xlim(0.2,8)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  xlab("v") +
  ylab("Log-likelihood") +
  ggtitle("v*=6")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white")
    # legend.position = "none"
  ) 
plot6
# 3*10
ggarrange(plot03, plot15, plot6, ncol = 3, nrow = 1) 

# tikz(file = "Selection nu.tex", standAlone=F,width = 7, height = 3)

# endoffile <- dev.off() 