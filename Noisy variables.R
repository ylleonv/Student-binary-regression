library(glmnet); library(ggplot2)

shape1  <- 2
shape2 <- 2

beta1 <- 1
alpha <- 0
amplitude_.5 <- 1
n_obs <- 250
n_test <- 250*0.2
n_training <- 250*0.8
p_var <- 5

s <-  0.05

n_sim <- 10

df_stu <- loglik_logit <- loglik_stu <- accu_logit <- accu_stu <- matrix(nrow = n_sim, ncol = 1)
out_coef_p <- out_coef_l <-  NULL

simu_f <- function(s,n_sim){
  # p_out <- 5
  # s <- 0.5
  # n_out <- n_training*p_out/100
  sim <- 0
  
  while (sim < n_sim) {
    # print(sim)
    sim <- sim+1
    z_beta <- rbeta(n = n_obs, shape1, shape2)
    Sigma <- matrix(0, nrow = p_var, ncol = p_var)
    diag(Sigma) <- 1
    Ex1 <- cbind(MASS::mvrnorm(n_obs, rep(0, p_var), Sigma))
    Ex <- (amplitude_.5-alpha)*((2*z_beta) - 1) / beta1[1]
    # Ex <- ((2*amplitude_.5*z_beta) / (beta1*(max(z_beta)-min(z_beta)))) + ((-amplitude_.5-alpha) / beta1)
    Ex <- cbind(Ex,Ex1[,-1])
    
    
    beta1 <- c(1, 0, 0, rep(0, p_var - 3))
    eta1 <- alpha + as.matrix(Ex) %*% (beta1)
    
    
    h_x_beta <- S(eta1,s)
    y1 <- rbinom(n = n_obs,size = 1,prob = h_x_beta)    
    
    data_sim <- as.data.frame(cbind(Ex, response = y1))
    dat_training <- data_sim[1:n_training,]
    dat_test <- data_sim[(n_training+1):n_obs,]
    
    # if (n_out > 0) {
    # n_out <- 2
    # dat_training[
    #   rownames(dat_training[dat_training$Ex %in% sort(dat_training[,"Ex"], F)[1:n_out],])
    #   ,"response"] <- rep(1,n_out)
    # # }
    
    mod_stu <- glm(response~.-1,family = binomial(Gosset(0.8)),data = dat_training, control = glm.control(maxit = 100))
    
    xfactors <- model.matrix(response~.-1, data = dat_training)
    lasso_stu <- cv.glmnet(xfactors, factor(dat_training$response), family = binomial(Gosset(0.8)), alpha = 1, lambda = NULL, intercept = F)
    lasso_stu_v <- data.frame(t(coef(lasso_stu, s = "lambda.min")[-1]))
    colnames(lasso_stu_v) <- rownames(coef(lasso_stu, s = "lambda.min"))[-1]
    
    lasso_log <- cv.glmnet(xfactors, factor(dat_training$response), family = binomial("logit"), alpha = 1, lambda = NULL, intercept = F)
    lasso_log_v <- data.frame(t(coef(lasso_log, s = "lambda.min")[-1]))
    colnames(lasso_log_v) <- rownames(coef(lasso_log, s = "lambda.min"))[-1]
    
    res1_l <- dplyr::bind_rows(((lasso_stu_v)),((lasso_log_v)))
    res1_l$link <- c("student","logit")  
    out_coef_l <- dplyr::bind_rows(out_coef_l,res1_l)
    
    loglik_stu[sim] <- logLik(mod_stu)
    step_stu <- step(mod_stu, direction = "both")
    p_stu <- data.frame(t(summary(step_stu)$coefficients[,4]))
    colnames(p_stu) <- rownames(summary(step_stu)$coefficients)
    
    mod_logit <- glm(response~.-1,family = binomial("logit"),data = dat_training, control = glm.control(maxit = 100))
    loglik_logit[sim] <- logLik(mod_logit)
    step_logit <- step(mod_logit, direction = "both")
    p_logit <- data.frame(t(summary(step_logit)$coefficients[,4]))
    colnames(p_logit) <- rownames(summary(step_logit)$coefficients)
    res1 <- dplyr::bind_rows(((p_stu)),((p_logit)))
    res1$link <- c("student","logit")    
    
    out_coef_p <- dplyr::bind_rows(out_coef_p,res1)
  }
  
  log_lik_res <- rbind(data.frame(loglik=loglik_logit,link="logistic",pi=s),
                       data.frame(loglik=loglik_stu,link="student",pi=s))
  
  return(list(log_lik_res,out_coef_p,out_coef_l))
}

p_s <- simu_f(0.05, 100)

log_lik_p_s <- p_s[[1]]
coef_p <- p_s[[2]]
coef_l <- p_s[[3]]

head(coef_p)
colnames(coef_p)[1] <- "V00"
library(ggplot2)

long_data <- pivot_longer(coef_p,cols = starts_with("V"))
long_data$link <- as.factor(long_data$link)
long_data$name <- as.factor(long_data$name)
head(long_data)

long_data %>% 
  group_by(link, name) %>% 
  dplyr::summarise(avg = mean(value, na.rm = T), cnt = sum(!is.na(value)), cnt_sig = sum(na.omit(value)<0.05) )

ggplot(data = long_data, mapping = aes(color = link, x = factor(name), y = value)) +
  # geom_violin(scale = "count")+
  geom_boxplot()+
  # geom_point()+
  geom_jitter(position=position_jitter(height=0, width=0.2))
  

head(coef_l)
colnames(coef_l)[1] <- "V00"
library(ggplot2)

long_data <- pivot_longer(coef_l,cols = starts_with("V"))
long_data$link <- as.factor(long_data$link)
long_data$name <- as.factor(long_data$name)
head(long_data)

long_data %>% 
  group_by(link, name) %>% 
  dplyr::summarise(avg = mean(value, na.rm = T), cnt = sum(value != 0) )

ggplot(data = long_data, mapping = aes(color = link, x = factor(name), y = value)) +
  # geom_violin(scale = "count")+
  geom_boxplot()+
  # geom_point()+
  geom_jitter(position=position_jitter(height=0, width=0.2))
  
  df_lik_p_s <- p_s[[3]]