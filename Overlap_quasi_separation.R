library(tidyr);library(dplyr); library(ggplot2); library(tikzDevice)

set.seed(5)

x1 <- rnorm(45,1,1)
x2 <- rnorm(45,0,1)

dat1 <- tbl_df(cbind(x1,x2))
# dat1[12,] <- t(c(1,0))
# dat1[13,] <- t(c(2,0.8))

vec1 <- c(25,2,5,11,21,23,27,28,32,33,34,35,36,41,44,6,8,17,42,19,40)

dat1$group <- 0
dat1[vec1,]$group <- 1

quasisep <- ggplot(dat1, aes(x=x1,y=x2,color=factor(group)))+
  # geom_text(
  #   label=rownames(dat1),
  #   check_overlap = T
  # )+
  geom_abline(aes(intercept = 0.92, slope = -0.824),linetype="dashed", colour = "gray")+
  geom_point()+
  scale_color_manual("", labels = c("0", "1"), values = c("#B4464B", "#B4AF46"))+
  # scale_shape_manual(values = c(1, 16))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none")
quasisep

head(dat1)
dat1$group <- as.factor(dat1$group)
library(GLMcat)
# glmcat(group~x1+x2, data = dat1, ratio = "reference",cdf = "logistic")
mod1 <- glm(group~x1+x2, data = dat1, family = binomial(), control = glm.control(trace = T))
summary(mod1)
mod1$linear.predictors
mod1$fitted.values
plot(mod1$linear.predictors,mod1$fitted.values)

separation <- ggplot(dat1[-c(11,12,31,10,2),], aes(x=x1,y=x2,color=factor(group)))+
  # geom_text(
  #   label=rownames(dat1),
  #   check_overlap = T
  # )+
  geom_abline(aes(intercept = 0.92, slope = -0.824),linetype="dashed", colour = "gray")+
  geom_point()+
  scale_color_manual("", labels = c("0", "1"), values = c("#B4464B", "#B4AF46"))+
  # scale_shape_manual(values = c(1, 16))+
  # geom_abline(aes(intercept = 1.2, slope = -1.05),linetype="dashed", colour = "gray")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none")
separation

mod1 <- glm(group~x1+x2, data = dat1[-c(11,12,31,10,2),], family = binomial(), control = glm.control(trace = T))
summary(mod1)
mod1$linear.predictors
mod1$fitted.values
plot(mod1$linear.predictors,mod1$fitted.values)

vec1 <- c(2,5,11,21,23,27,28,32,33,35,36,41,44,8,17,42,19,40,15,4)

dat1$group <- 0
dat1[vec1,]$group <- 1

overlap <- ggplot(dat1, aes(x=x1,y=x2,color=factor(group)))+
  # geom_text(
  #   label=rownames(dat1),
  #   check_overlap = T
  # )+
  geom_point()+
  scale_color_manual("", labels = c("0", "1"), values = c("#B4464B", "#B4AF46"))+
  # scale_shape_manual(values = c(1, 16))+
  # geom_abline(aes(intercept = 0.92, slope = -0.824),linetype="dashed", colour = "gray")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none")
overlap

mod1 <- glm(group~x1+x2, data = dat1, family = binomial(), control = glm.control(trace = T))
summary(mod1)
mod1$linear.predictors
mod1$fitted.values
plot(mod1$linear.predictors,mod1$fitted.values)

# PDF 8*3
# tikz(file = "Overlap_quasi_separation.tex", standAlone=F, width = 6.7, height = 2.5)
# endoffile <- dev.off() 

library(ggpubr)
ggarrange(separation,quasisep,overlap,ncol = 3, nrow = 1, labels = c("",""))
# 4*10
# library(shiny)
# library(ggplot2)
# library(Cairo)   # For nicer ggplot2 output when deployed on Linux
# 
# ui <- fluidPage(
#   fluidRow(
#     column(width = 4,
#            plotOutput("plot1", height = 300,
#                       # Equivalent to: click = clickOpts(id = "plot_click")
#                       click = "plot1_click",
#                       brush = brushOpts(
#                         id = "plot1_brush"
#                       )
#            )
#     )
#   ),
#   fluidRow(
#     column(width = 6,
#            h4("Points near click"),
#            verbatimTextOutput("click_info")
#     ),
#     column(width = 6,
#            h4("Brushed points"),
#            verbatimTextOutput("brush_info")
#     )
#   )
# )
# 
# server <- function(input, output) {
#   output$plot1 <- renderPlot({
#     set.seed(5)
#     
#     x1 <- rnorm(45,1,1)
#     x2 <- rnorm(45,0,1)
#     
#     dat1 <- data.frame(cbind(x1,x2))
#     # dat1[12,] <- t(c(1,0))
#     # dat1[13,] <- t(c(2,0.8))
#     
#     ggplot(dat1, aes(x=x1,y=x2))+
#       geom_text(
#         label=rownames(dat1),
#         check_overlap = T
#       )+
#       geom_abline(aes(intercept = 0.9, slope = -0.835))
#   })
#   
#   output$click_info <- renderPrint({
#     # Because it's a ggplot2, we don't need to supply xvar or yvar; if this
#     # were a base graphics plot, we'd need those.
#     nearPoints(dat1, input$plot1_click, addDist = TRUE)
#   })
#   
#   output$brush_info <- renderPrint({
#     brushedPoints(dat1, input$plot1_brush)
#   })
# }
# 
# shinyApp(ui, server)
# 
# 
# library(GLMcat)
# data(TravelChoice)
# debug(glmcat)
# discrete_cm(formula = choice ~ hinc[car] + gc + invt,
#             case_id = "indv",alternatives = "mode", reference = "air",
#             data = TravelChoice, alternative_specific = c("gc", "invt"),
#             cdf = "logistic")
