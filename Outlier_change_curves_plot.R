new_vaso2 <- rbind(vaso, c(-1,0,1), c(-3,0,1), c(-5,0,1))

mod_01 <- glm(factor(vaso)~vol,  
             data = new_vaso2[-c(40:42),],
             family = binomial(link = Gosset(0.6)))
mod_02 <- glm(factor(vaso)~vol,  
             data = new_vaso2[-c(40:42),],
             family = binomial(link = "logit"))

mod_1 <- glm(factor(vaso)~vol,  
             data = new_vaso2[-c(41:42),],
             family = binomial(link = Gosset(0.6)))
mod_2 <- glm(factor(vaso)~vol,  
             data = new_vaso2[-c(41:42),],
             family = binomial(link = "logit"))

mod_3 <- glm(factor(vaso)~vol,  
             data = new_vaso2[-c(40,42),],
             family = binomial(link = Gosset(0.6)))
mod_4 <- glm(factor(vaso)~vol,  
             data = new_vaso2[-c(40,42),],
             family = binomial(link = "logit"))

mod_5 <- glm(factor(vaso)~vol,  
             data = new_vaso2[-c(40,41),],
             family = binomial(link = Gosset(0.6)))
mod_6 <- glm(factor(vaso)~vol,  
             data = new_vaso2[-c(40,41),],
             family = binomial(link = "logit"))


dat1 <- data.frame(vol = seq(-5,1.5,by=0.1))


new_vaso2$new <- "o"
new_vaso2$new[40] <- "n1"
new_vaso2$new[41] <- "n2"
new_vaso2$new[42] <- "n3"

plot_logi <- ggplot(newdata = new_vaso2) +
  geom_text(data=new_vaso2, aes(x = vol, y = vaso, label=vaso, color = new),
            check_overlap = T)+
  # geom_point(data=new_vaso2, aes(x = vol, y = vaso, label=vaso, color = new),
            # check_overlap = T)+
  scale_colour_manual(values=c("#FF5733", "#52854C", "#4E84C4", "#000000"))+
  # geom_hline(yintercept = 0.5, color = "gray", linetype = "longdash", size = 1)+
  stat_function(fun = function(x) {0.5}, color = "gray", linetype = "dotted", size = 1)+
  # geom_path(aes(y=as.vector(predict(mod_1, newdata = data.frame(vol=dat1), type = "response")),
  #           x = as.vector(dat1$vol)),linetype = "longdash", color = "#FF5733")+
  geom_path(aes(y=as.vector(predict(mod_2, newdata = data.frame(vol=dat1), type = "response")),
                x = as.vector(dat1$vol)),linetype = "solid", color = "#FF5733")+
  # geom_path(aes(y=as.vector(predict(mod_3, newdata = data.frame(vol=dat1), type = "response")),
  #               x = as.vector(dat1$vol)),linetype = "longdash", color = "#52854C")+
  geom_path(aes(y=as.vector(predict(mod_4, newdata = data.frame(vol=dat1), type = "response")),
                x = as.vector(dat1$vol)),linetype = "solid", color = "#52854C")+
  # geom_path(aes(y=as.vector(predict(mod_5, newdata = data.frame(vol=dat1), type = "response")),
  #               x = as.vector(dat1$vol)),linetype = "longdash", color = "#4E84C4")+
  geom_path(aes(y=as.vector(predict(mod_6, newdata = data.frame(vol=dat1), type = "response")),
                x = as.vector(dat1$vol)),linetype = "solid", color = "#4E84C4")+
  # geom_path(aes(y=as.vector(predict(mod_01, newdata = data.frame(vol=dat1), type = "response")),
  #               x = as.vector(dat1$vol)),linetype = "longdash", color = "#000000")+
  geom_path(aes(y=as.vector(predict(mod_02, newdata = data.frame(vol=dat1), type = "response")),
                x = as.vector(dat1$vol)),linetype = "solid", color = "#000000")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position = "none") +
  scale_y_continuous("", c(0,0.5,1), labels=c("0","0.5","1")) +
  xlab("ln of Volume") + ylab("") + ggtitle("Logit link")

plot_stu <- ggplot(newdata = new_vaso2) +
  geom_text(data=new_vaso2, aes(x = vol, y = vaso, label=vaso, color = new),
            check_overlap = T)+
  # geom_point(data=new_vaso2, aes(x = vol, y = vaso, label=vaso, color = new),
  # check_overlap = T)+
  scale_colour_manual(values=c("#FF5733", "#52854C", "#4E84C4", "#000000"))+
  # geom_hline(yintercept = 0.5, color = "gray", linetype = "longdash", size = 1)+
  stat_function(fun = function(x) {0.5}, color = "gray", linetype = "dotted", size = 1)+
  geom_path(aes(y=as.vector(predict(mod_1, newdata = data.frame(vol=dat1), type = "response")),
            x = as.vector(dat1$vol)),linetype = "longdash", color = "#FF5733")+
  # geom_path(aes(y=as.vector(predict(mod_2, newdata = data.frame(vol=dat1), type = "response")),
  #               x = as.vector(dat1$vol)),linetype = "solid", color = "#FF5733")+
  geom_path(aes(y=as.vector(predict(mod_3, newdata = data.frame(vol=dat1), type = "response")),
                x = as.vector(dat1$vol)),linetype = "longdash", color = "#52854C")+
  # geom_path(aes(y=as.vector(predict(mod_4, newdata = data.frame(vol=dat1), type = "response")),
  #               x = as.vector(dat1$vol)),linetype = "solid", color = "#52854C")+
  geom_path(aes(y=as.vector(predict(mod_5, newdata = data.frame(vol=dat1), type = "response")),
                x = as.vector(dat1$vol)),linetype = "longdash", color = "#4E84C4")+
  # geom_path(aes(y=as.vector(predict(mod_6, newdata = data.frame(vol=dat1), type = "response")),
  #               x = as.vector(dat1$vol)),linetype = "solid", color = "#4E84C4")+
  geom_path(aes(y=as.vector(predict(mod_01, newdata = data.frame(vol=dat1), type = "response")),
                x = as.vector(dat1$vol)),linetype = "longdash", color = "#000000")+
  # geom_path(aes(y=as.vector(predict(mod_02, newdata = data.frame(vol=dat1), type = "response")),
  #               x = as.vector(dat1$vol)),linetype = "solid", color = "#000000")
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position = "none") +
  scale_y_continuous("", c(0,0.5,1), labels=c("0","0.5","1")) +
  xlab("ln of Volume") + ylab("") + ggtitle("Student link")

ggarrange(plot_logi, plot_stu, ncol = 2, nrow = 1) 

# 4*11
# tikz(file = "Outlier_change_curves_plot.tex", standAlone=F,width = 7, height = 3)
# endoffile <- dev.off() 



