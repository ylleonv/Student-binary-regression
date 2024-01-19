shape1  <- 2
shape2 <- 2

beta1 <- 1
alpha <- 0
# amplitude_.5 <- 5
n_obs <- 100
n_test <- 100*0.2
n_training <- 100*0.8
p_var <- 1

z_beta <- rbeta(n = n_obs, shape1, shape2)
z_beta <- rnorm(n_obs,0,1)
Ex <- z_beta
# Sigma <- matrix(0.5, nrow = p_var, ncol = p_var)
# diag(Sigma) <- 1
# Ex1 <- cbind(mvrnorm(n_obs, rep(0, p_var), Sigma))
# 
# Ex <- (amplitude_.5-alpha)*((2*z_beta) - 1) / beta1
# summary(Ex)
# Ex <- cbind(Ex,Ex1[,-1])

beta1 <- c(1, rep(0, p_var - 1))

eta1 <- alpha + as.matrix(Ex) %*% (beta1)

S<-function(x,s){ return  (1/(1+exp(-x/s)))}

s <- c(0.05,0.1,0.3,0.8)

h_x_beta <- S(eta1,s[1])
y1 <- rbinom(n = n_obs,size = 1,prob = h_x_beta)    
data_sim1 <- as.data.frame(cbind(Ex, response = y1))
data_sim1$New <- s[1]

for (variable in 2:length(s)) {
  h_x_beta <- S(eta1,s[variable])
  y1 <- rbinom(n = n_obs,size = 1,prob = h_x_beta)    
  data_sim <- as.data.frame(cbind(Ex, response = y1))
  data_sim$New <- s[variable]
  data_sim1 <- rbind(data_sim1,data_sim)
}

summary(data_sim1)

mean_df <- data_sim1 %>% 
  group_by(New, response) %>% 
  summarise(min = min(Ex))

mean_df <- mean_df[mean_df$response == 1,]
mean_df <- as.data.frame(mean_df)
summary(mean_df)

max_df <- data_sim1 %>% 
  group_by(New, response) %>% 
  summarise(max = max(Ex))

max_df <- max_df[max_df$response == 0,]
max_df <- as.data.frame(max_df)
summary(max_df)

ggplot(data = data_sim1)+
  geom_boxplot( aes(x=Ex, fill = factor(response))) +
  # geom_vline(xintercept = min(Ex), linetype="dashed")+
  geom_vline(data = mean_df, aes(xintercept = min), linetype="dashed", color = "#B4AF46")+
  geom_vline(data = max_df, aes(xintercept = max), linetype="dashed", color = "#B4464B")+
  facet_grid(.~New)+
  theme(
    strip.background = element_rect(
      fill="white"
    ),
    # panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  ) +
  # scale_color_manual("", labels = c("0", "1"), values = c("#4E84C4", "#52854C"))+
  # scale_fill_manual("", labels = c("0", "1"), values = c("#4E84C4", "#52854C"))
  scale_fill_manual("", labels = c("0", "1"), values = c("#B4464B", "#B4AF46"))+
  scale_x_continuous(name="x", breaks=seq(-5, 5, 5))+
  scale_y_continuous(labels = NULL, breaks = NULL) 

# pdf 12*3

# PDF 8*3
# tikz(file = "Plot separation.tex", standAlone=F, width = 7, height = 2.5)
# endoffile <- dev.off() 

head(data_sim1)
dat_out <- data_sim1[data_sim1$New==0.3,]
dat_out$Out <- "no"
dat_out[dat_out$Ex==min(dat_out$Ex),"response"] <- 1
dat_out[dat_out$Ex==min(dat_out$Ex),"Out"] <- "yes"
str(dat_out)
dat_out[nrow(dat_out)+1,] <- c(-2.8,1,0.3,"yes")
dat_out$Ex <- as.numeric(dat_out$Ex)
dat_out$response <- as.numeric(dat_out$response)
dat_out$New <- as.numeric(dat_out$New)

ggplot(data = dat_out)+
  geom_boxplot( aes(x=Ex,y=response, color = factor(response)), outlier.shape = NA) +
  geom_jitter( aes(x=Ex,y=response, color = factor(response), shape = factor(Out)), height = 0.2, width = 0) +
  # geom_vline(xintercept = min(Ex), linetype="dashed")+
  # geom_vline(data = mean_df, aes(xintercept = min), linetype="dashed", color = "#B4AF46")+
  # geom_vline(data = max_df, aes(xintercept = max), linetype="dashed", color = "#B4464B")+
  theme(
    strip.background = element_rect(
      fill="white"
    ),
    # panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "white", colour = "grey50"),
    legend.key = element_rect(fill = "white"),
    legend.position = "none"
  ) +
  # scale_color_manual("", labels = c("0", "1"), values = c("#4E84C4", "#52854C"))+
  scale_shape_manual("", labels = c("0", "1"), values = c(1, 8))+
  scale_color_manual("", labels = c("0", "1"), values = c("#B4464B", "#B4AF46"))+
  scale_x_continuous(name="$x$", breaks=seq(-5, 5, 5))+
  scale_y_continuous(name="$y$", breaks=c(0,1)) 

tikz(file = "tex/Plot_separation_outliers.tex", standAlone=F, width = 5, height = 3)
endoffile <- dev.off()

