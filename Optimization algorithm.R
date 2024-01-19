# Functions optimization --------------------------------------------------

g_deviance <- function(nu, data_in) {
  print(nu)
  logLik(glm(y~.,family = binomial(Gosset(nu)),data = data_in, maxit = 100))
}

S<-function(x,s){ return  (1/(1+exp(-x/s)))}

n_obs <- 100
s <- 0.1
n_obs <- n_obs*1.3

alpha <- 0
beta1 <- 1
Sigma <- matrix(0, nrow = 1, ncol = 1)
diag(Sigma) <- 1
Ex <- cbind(MASS::mvrnorm(n_obs, rep(0, 1), Sigma))
eta1 <- alpha + as.matrix(Ex) %*% (beta1)
h_x_beta <- S(eta1,s)
# h_x_beta <- pt(eta1,s)
y <- rbinom(n = n_obs,size = 1,prob = h_x_beta)  
dat2 <- data.frame(y,Ex)
dat2$y <- as.factor(dat2$y)
dat3 <- dat2
n_out_f <- 2
outliers <- data.frame(y=rep(0,n_out_f),data.frame(MASS::mvrnorm(n = n_out_f, rep(2, 1), Sigma))) 
colnames(outliers) <- colnames(dat2)
dat3 <- rbind(dat2,outliers)
head(dat3)

library(ggplot2)
ggplot(data = dat3,mapping = aes(x=Ex,y=y))+
  geom_point()

a2 <- optimize(g_deviance, c(0.25, 1), maximum = T, dat3, tol = 0.1)
# undebug(optimize)
# a2$maximum
# a2$objective
# 
# lower <- 0.25
# upper <- 1
# 
# phi <- ((sqrt(5) - 1)/2)
# x_1 = lower + (1-phi)*(upper-lower) 
# 
# x_1
# 
# x_2 = lower + phi*(upper-lower)
# 
# 
# logLik(glm(y~.,family = binomial(Gosset(x_1)),data = dat3, maxit = 100))
# logLik(glm(y~.,family = binomial(Gosset(x_2)),data = dat3, maxit = 100))
# 
# 
# lower <- 0.25
# upper <- x_1
# 
# x_1 = lower + (1-phi)*(upper-lower) 
# x_1
# x_2 = lower + phi*(upper-lower)
# 
# logLik(glm(y~.,family = binomial(Gosset(x_1)),data = dat3, maxit = 100))
# logLik(glm(y~.,family = binomial(Gosset(x_2)),data = dat3, maxit = 100))
# 
# lower <- x_1
# upper <- 0.5364745
# 
# x_1 = lower + (1-phi)*(upper-lower) 
# x_1
# x_2 = lower + phi*(upper-lower)
# x_2
