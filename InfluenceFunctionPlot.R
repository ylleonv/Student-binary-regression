
# Influence function plot -------------------------------------------------

eta <- seq(-6,6,0.01)
logis <- dlogis(eta)*eta/(plogis(eta)*(1-plogis(eta)))
normal <- dnorm(eta)*eta/(pnorm(eta)*(1-pnorm(eta)))
cauchy <- dcauchy(eta)*eta/(pcauchy(eta)*(1-pcauchy(eta)))
student_1 <- dt(eta,df = 1)*eta/(pt(eta,df = 1)*(1-pt(eta,df = 1)))
student_03 <- dt(eta,df = 0.3)*eta/(pt(eta,df = 0.3)*(1-pt(eta,df = 0.3)))
# student_08 <- dt(eta,df = 0.8)*eta/(pt(eta,df = 0.8)*(1-pt(eta,df = 0.8)))
student_2 <- dt(eta,df = 2)*eta/(pt(eta,df = 2)*(1-pt(eta,df = 2)))
student_4 <- dt(eta,df = 4)*eta/(pt(eta,df = 4)*(1-pt(eta,df = 4)))
# gompertz <- (exp(eta-exp(eta)))*eta/((1-exp(-exp(eta)))*(1-(1-exp(-exp(eta)))))
# gumbel <- (exp(-(eta+exp(-eta))))*eta/((exp(-exp(-eta)))*(1-(exp(-exp(-eta)))))


links <- rbind(data.frame(eta,value=logis,Link="Logistic", Dis1 = "Logistic"),
               data.frame(eta,value=normal,Link="Normal", Dis1 = "Normal"),
               data.frame(eta,value=student_1,Link="Student", Dis1 = "Student(1)"),
               data.frame(eta,value=student_03,Link="Student", Dis1 = "Student(0.3)"),
               # data.frame(eta,value=student_08,Link="Student", Dis1 = "Student(0.8)"),
               # data.frame(eta,value=cauchy,Link="Cauchy", Dis1 = "Cauchy"),
               # data.frame(eta,value=gumbel,Link="Gumbel", Dis1 = "Gumbel"),
               data.frame(eta,value=student_2,Link="Student", Dis1 = "Student(2)"),
               data.frame(eta,value=student_4,Link="Student", Dis1 = "Student(4)"))

summary(links)

library(RColorBrewer)

ggplot(data = links,aes(x = eta, y = value, color=factor(Dis1), linetype = factor(Dis1)))+
  # geom_line(aes(x = eta, y = value, color=factor(Link), linetype = factor(Link))) +
  # labs(title="", x="eta", y="s(eta)", color = "CDF")+
  ylim(-6,6)+ylab("")+xlim(-6,6)+
  geom_line(lwd = 1) +
  scale_linetype_manual(values = c(rep("solid", 2), "twodash","longdash", "dotdash", "dashed")) +
  scale_color_manual(values = c("#E69F00", "#C70039","#00AFBB","#00AFBB","#00AFBB","#00AFBB"))+
  # scale_color_discrete(name="") +
  theme(  strip.background = element_rect(
      fill="yellow"
    ),
    # panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "yellow", colour = "grey50"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "yellow", colour = "yellow"),
    legend.background = element_rect(fill = "yellow"),
    legend.key = element_rect(fill = NA)
        )

# save 6*8
tikz(file = "INFLUENCE_FUNCTION.tex", standAlone=F,width = 3.5, height = 2)
endoffile <- dev.off()

