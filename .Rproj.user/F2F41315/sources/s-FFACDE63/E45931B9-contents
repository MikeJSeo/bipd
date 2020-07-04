# Loading
library("ggplot2")
library("readxl")
library("grid")
require("gridExtra")

#setwd("C:/Users/mike/Desktop/Github/phd/varselect/simulation_results_renumbered")
setwd("~/GitHub/phd/varselect/simulation_results")

simulation1 <- read_excel("simulation1.result.xlsx")
simulation2 <- read_excel("simulation2.result.xlsx")
simulation3 <- read_excel("simulation3.result.xlsx")
simulation4 <- read_excel("simulation4.result.xlsx")
simulation5 <- read_excel("simulation5.result.xlsx")
simulation6 <- read_excel("simulation6.result.xlsx")
simulation7 <- read_excel("simulation7.result.xlsx")
simulation8 <- read_excel("simulation8.result.xlsx")
simulation9 <- read_excel("simulation9.result.xlsx")
simulation10 <- read_excel("simulation10.result.xlsx")
simulation11 <- read_excel("simulation11.result.xlsx")
simulation12 <- read_excel("simulation12.result.xlsx")
simulation13 <- read_excel("simulation13.result.xlsx")
simulation14 <- read_excel("simulation14.result.xlsx")
simulation15 <- read_excel("simulation15.result.xlsx")
simulation16 <- read_excel("simulation16.result.xlsx")
simulation17 <- read_excel("simulation17.result.xlsx")
simulation18 <- read_excel("simulation18.result.xlsx")
simulation19 <- read_excel("simulation19.result.xlsx")
simulation20 <- read_excel("simulation20.result.xlsx")
simulation21 <- read_excel("simulation21.result.xlsx")
simulation22 <- read_excel("simulation22.result.xlsx")
simulation23 <- read_excel("simulation23.result.xlsx")
simulation24 <- read_excel("simulation24.result.xlsx")
simulation25 <- read_excel("simulation25.result.xlsx")
simulation26 <- read_excel("simulation26.result.xlsx")
simulation27 <- read_excel("simulation27.result.xlsx")
simulation28 <- read_excel("simulation28.result.xlsx")
simulation29 <- read_excel("simulation29.result.xlsx")
simulation30 <- read_excel("simulation30.result.xlsx")
simulation31 <- read_excel("simulation31.result.xlsx")
simulation32 <- read_excel("simulation32.result.xlsx")
simulation33 <- read_excel("simulation33.result.xlsx")
simulation34 <- read_excel("simulation34.result.xlsx")
simulation35 <- read_excel("simulation35.result.xlsx")
simulation36 <- read_excel("simulation36.result.xlsx")
simulation37 <- read_excel("simulation37.result.xlsx")
simulation38 <- read_excel("simulation38.result.xlsx")
simulation39 <- read_excel("simulation39.result.xlsx")
simulation40 <- read_excel("simulation40.result.xlsx")
simulation41 <- read_excel("simulation41.result.xlsx")
simulation42 <- read_excel("simulation42.result.xlsx")
simulation43 <- read_excel("simulation43.result.xlsx")
simulation44 <- read_excel("simulation44.result.xlsx")
simulation45 <- read_excel("simulation45.result.xlsx")
simulation46 <- read_excel("simulation46.result.xlsx")
simulation47 <- read_excel("simulation47.result.xlsx")
simulation48 <- read_excel("simulation48.result.xlsx")
simulation49 <- read_excel("simulation49.result.xlsx")
simulation50 <- read_excel("simulation50.result.xlsx")
simulation51 <- read_excel("simulation51.result.xlsx")
simulation52 <- read_excel("simulation52.result.xlsx")
simulation53 <- read_excel("simulation53.result.xlsx")
simulation54 <- read_excel("simulation54.result.xlsx")
simulation55 <- read_excel("simulation55.result.xlsx")
simulation56 <- read_excel("simulation56.result.xlsx")
simulation57 <- read_excel("simulation57.result.xlsx")
simulation58 <- read_excel("simulation58.result.xlsx")
simulation59 <- read_excel("simulation59.result.xlsx")
simulation60 <- read_excel("simulation60.result.xlsx")


make_data <- function(simulation_a, simulation_b, a_name, b_name, mse = 4, xlab){
  
  
  data1 <- as.data.frame(simulation_a)[,mse, drop = FALSE]
  colnames(data1) <- "error"
  data1$models <- c("A", "B", "C", "D", "E", "F", "G")
  data1$simulations <- a_name
  
  data2 <- as.data.frame(simulation_b)[,mse, drop = FALSE]
  colnames(data2) <- "error"
  data2$models <- c("A", "B", "C", "D", "E", "F", "G")
  data2$simulations <- b_name
  
  data <- rbind(data1, data2)
  
  data$simulations <- factor(c(rep(a_name, 7), rep(b_name,7)), levels = c(a_name, b_name))
  list(data = data, xlab = xlab)
}


data1 <- make_data(simulation1, simulation16, a_name = "Scenario 1 \n N = 5", b_name = "Scenario 16 \n N = 10",
                  xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.2, no effect modification")

data2 <- make_data(simulation2, simulation17, a_name = "Scenario 2 \n N = 5", b_name = "Scenario 17 \n N = 10",
                   xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.5, no effect modification")

data3 <- make_data(simulation3, simulation18, a_name = "Scenario 3 \n N = 5", b_name = "Scenario 18 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, small effect modification")

data4 <- make_data(simulation4, simulation19, a_name = "Scenario 4 \n N = 5", b_name = "Scenario 19 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, small effect modification")

data5 <- make_data(simulation5, simulation20, a_name = "Scenario 5 \n N = 5", b_name = "Scenario 20 \n N = 10", 
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, large effect modification")

data6 <- make_data(simulation6, simulation21, a_name = "Scenario 6 \n N = 5", b_name = "Scenario 21 \n N = 10", 
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, large effect modification")

data7 <- make_data(simulation7, simulation22, a_name = "Scenario 7 \n N = 5", b_name = "Scenario 22 \n N = 10",
                   xlab = "10 covariates, 10 effect modifiers \n \u03c4 = 0.2, small effect modification")

data8 <- make_data(simulation8, simulation23, a_name = "Scenario 8 \n N = 5", b_name = "Scenario 23 \n N = 10", 
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, small effect modification")

data9 <- make_data(simulation9, simulation24, a_name = "Scenario 9 \n N = 5", b_name = "Scenario 24 \n N = 10", 
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, small effect modification")

data10 <- make_data(simulation10, simulation25, a_name = "Scenario 10 \n N = 5", b_name = "Scenario 25 \n N = 10", 
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, large effect modification")

data11 <- make_data(simulation11, simulation26, a_name = "Scenario 11 \n N = 5", b_name = "Scenario 26 \n N = 10", 
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, large effect modification")

data12 <- make_data(simulation12, simulation27, a_name = "Scenario 12 \n N = 5", b_name = "Scenario 27 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, small effect modification")

data13 <- make_data(simulation13, simulation28, a_name = "Scenario 13 \n N = 5", b_name = "Scenario 28 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, small effect modification")

data14 <- make_data(simulation14, simulation29, a_name = "Scenario 14 \n N = 5", b_name = "Scenario 29 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, large effect modification")

data15 <- make_data(simulation15, simulation30, a_name = "Scenario 15 \n N = 5", b_name = "Scenario 30 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, large effect modification")

## binary

data16 <- make_data(simulation31, simulation46, a_name = "Scenario 31 \n N = 5", b_name = "Scenario 46 \n N = 10",
                   xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.2, no effect modification")

data17 <- make_data(simulation32, simulation47, a_name = "Scenario 32 \n N = 5", b_name = "Scenario 47 \n N = 10",
                   xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.5, no effect modification")

data18 <- make_data(simulation33, simulation48, a_name = "Scenario 33 \n N = 5", b_name = "Scenario 48 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, small effect modification")

data19 <- make_data(simulation34, simulation49, a_name = "Scenario 34 \n N = 5", b_name = "Scenario 49 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, small effect modification")

data20 <- make_data(simulation35, simulation50, a_name = "Scenario 35 \n N = 5", b_name = "Scenario 50 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, large effect modification")

data21 <- make_data(simulation36, simulation51, a_name = "Scenario 36 \n N = 5", b_name = "Scenario 51 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, large effect modification")

data22 <- make_data(simulation37, simulation52, a_name = "Scenario 37 \n N = 5", b_name = "Scenario 52 \n N = 10",
                   xlab = "10 covariates, 10 effect modifiers \n \u03c4 = 0.5, large effect modification")

data23 <- make_data(simulation38, simulation53, a_name = "Scenario 38 \n N = 5", b_name = "Scenario 53 \n N = 10",
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, small effect modification")

data24 <- make_data(simulation39, simulation54, a_name = "Scenario 39 \n N = 5", b_name = "Scenario 54 \n N = 10",
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, small effect modification")

data25 <- make_data(simulation40, simulation55, a_name = "Scenario 40 \n N = 5", b_name = "Scenario 55 \n N = 10",
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, large effect modification")

data26 <- make_data(simulation41, simulation56, a_name = "Scenario 41 \n N = 5", b_name = "Scenario 56 \n N = 10",
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, large effect modification")

data27 <- make_data(simulation42, simulation57, a_name = "Scenario 42 \n N = 5", b_name = "Scenario 57 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, small effect modification")

data28 <- make_data(simulation43, simulation58, a_name = "Scenario 43 \n N = 5", b_name = "Scenario 58 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, small effect modification")

data29 <- make_data(simulation44, simulation59, a_name = "Scenario 44 \n N = 5", b_name = "Scenario 59 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, large effect modification")

data30 <- make_data(simulation45, simulation60, a_name = "Scenario 45 \n N = 5", b_name = "Scenario 60 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, large effect modification")




plot1 <- ggplot(data=data1$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data1$xlab)

plot2 <- ggplot(data=data2$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data2$xlab)

plot3 <- ggplot(data=data3$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data3$xlab)

plot4 <- ggplot(data=data4$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data4$xlab)

plot5 <- ggplot(data=data5$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data5$xlab)

plot6 <- ggplot(data=data6$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data6$xlab)

plot7 <- ggplot(data=data7$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data7$xlab)

plot8 <- ggplot(data=data8$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data8$xlab)

plot9 <- ggplot(data=data9$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data9$xlab)

plot10 <- ggplot(data=data10$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data10$xlab)

plot11 <- ggplot(data=data11$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data11$xlab)

plot12 <- ggplot(data=data12$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data12$xlab)

plot13 <- ggplot(data=data13$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data13$xlab)

plot13 <- ggplot(data=data13$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data13$xlab)

plot14 <- ggplot(data=data14$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data14$xlab)

plot15 <- ggplot(data=data15$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data15$xlab)

plot16 <- ggplot(data=data16$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data16$xlab)

plot17 <- ggplot(data=data17$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data17$xlab)

plot18 <- ggplot(data=data18$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data18$xlab)

plot19 <- ggplot(data=data19$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data19$xlab)

plot20 <- ggplot(data=data20$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data20$xlab)

plot21 <- ggplot(data=data21$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data21$xlab)

plot22 <- ggplot(data=data22$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data22$xlab)

plot23 <- ggplot(data=data23$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data23$xlab)

plot24 <- ggplot(data=data24$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data24$xlab)

plot25 <- ggplot(data=data25$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data25$xlab)

plot26 <- ggplot(data=data26$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data26$xlab)

plot27 <- ggplot(data=data27$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data27$xlab)

plot28 <- ggplot(data=data28$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data28$xlab)

plot29 <- ggplot(data=data29$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data29$xlab)

plot30 <- ggplot(data=data30$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data30$xlab)



windowsFonts(Times = windowsFont("Times New Roman"))
fg <- frameGrob()
tg <- textGrob("Continuous outcome", gp = gpar(fontsize = 28, fontfamily= "Times"))
rg <- rectGrob(x = tg$x+unit(4, "mm"), y = tg$y, width = stringWidth(tg$label)*5 + unit(9, "mm") ,                 
               height = stringHeight(tg$label) + unit(10,"mm"), gp = gpar(fill = "light grey", lty = 0))
fg <- packGrob(fg, rg)
fg <- packGrob(fg, tg)


fg2 <- frameGrob()
tg2 <- textGrob("Binary outcome", gp = gpar(fontsize = 28, fontfamily= "Times"))
rg2 <- rectGrob(x = tg2$x+unit(4, "mm"), y = tg2$y, width = stringWidth(tg$label)*5 + unit(9, "mm"),                 
               height = stringHeight(tg2$label) + unit(10,"mm"), gp = gpar(fill = "light grey", lty = 0))
fg2 <- packGrob(fg2, rg2)
fg2 <- packGrob(fg2, tg2)

grid.arrange(arrangeGrob(plot1,plot2, plot3,plot4,plot5,plot6,plot7,plot8,plot9, plot10, plot11, plot12, plot13, plot14, plot15,
                         top=fg, ncol=3, as.table = FALSE),
             arrangeGrob(plot16,plot17,plot18,plot19,plot20,plot21,plot22,plot23,plot24,plot25,plot26,plot27,plot28,plot29,plot30,
                         top=fg2, ncol=3, as.table = FALSE),
             left = textGrob("Patient specific treatment MSE", rot = 90, vjust = 0.5), ncol=2)




################################ treatment MSE


data1 <- make_data(simulation1, simulation16, a_name = "Scenario 1 \n N = 5", b_name = "Scenario 16 \n N = 10",
                   xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.2, no effect modification", mse = 3)

data2 <- make_data(simulation2, simulation17, a_name = "Scenario 2 \n N = 5", b_name = "Scenario 17 \n N = 10",
                   xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.5, no effect modification", mse = 3)

data3 <- make_data(simulation3, simulation18, a_name = "Scenario 3 \n N = 5", b_name = "Scenario 18 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, small effect modification", mse = 3)

data4 <- make_data(simulation4, simulation19, a_name = "Scenario 4 \n N = 5", b_name = "Scenario 19 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, small effect modification", mse = 3)

data5 <- make_data(simulation5, simulation20, a_name = "Scenario 5 \n N = 5", b_name = "Scenario 20 \n N = 10", 
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, large effect modification", mse = 3)

data6 <- make_data(simulation6, simulation21, a_name = "Scenario 6 \n N = 5", b_name = "Scenario 21 \n N = 10", 
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, large effect modification", mse = 3)

data7 <- make_data(simulation7, simulation22, a_name = "Scenario 7 \n N = 5", b_name = "Scenario 22 \n N = 10",
                   xlab = "10 covariates, 10 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 3)

data8 <- make_data(simulation8, simulation23, a_name = "Scenario 8 \n N = 5", b_name = "Scenario 23 \n N = 10", 
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 3)

data9 <- make_data(simulation9, simulation24, a_name = "Scenario 9 \n N = 5", b_name = "Scenario 24 \n N = 10", 
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, small effect modification", mse = 3)

data10 <- make_data(simulation10, simulation25, a_name = "Scenario 10 \n N = 5", b_name = "Scenario 25 \n N = 10", 
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, large effect modification", mse = 3)

data11 <- make_data(simulation11, simulation26, a_name = "Scenario 11 \n N = 5", b_name = "Scenario 26 \n N = 10", 
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 3)

data12 <- make_data(simulation12, simulation27, a_name = "Scenario 12 \n N = 5", b_name = "Scenario 27 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 3)

data13 <- make_data(simulation13, simulation28, a_name = "Scenario 13 \n N = 5", b_name = "Scenario 28 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, small effect modification", mse = 3)

data14 <- make_data(simulation14, simulation29, a_name = "Scenario 14 \n N = 5", b_name = "Scenario 29 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, large effect modification", mse = 3)

data15 <- make_data(simulation15, simulation30, a_name = "Scenario 15 \n N = 5", b_name = "Scenario 30 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 3)

## binary

data16 <- make_data(simulation31, simulation46, a_name = "Scenario 31 \n N = 5", b_name = "Scenario 46 \n N = 10",
                    xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.2, no effect modification", mse = 3)

data17 <- make_data(simulation32, simulation47, a_name = "Scenario 32 \n N = 5", b_name = "Scenario 47 \n N = 10",
                    xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.5, no effect modification", mse = 3)

data18 <- make_data(simulation33, simulation48, a_name = "Scenario 33 \n N = 5", b_name = "Scenario 48 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, small effect modification", mse = 3)

data19 <- make_data(simulation34, simulation49, a_name = "Scenario 34 \n N = 5", b_name = "Scenario 49 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, small effect modification", mse = 3)

data20 <- make_data(simulation35, simulation50, a_name = "Scenario 35 \n N = 5", b_name = "Scenario 50 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, large effect modification", mse = 3)

data21 <- make_data(simulation36, simulation51, a_name = "Scenario 36 \n N = 5", b_name = "Scenario 51 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, large effect modification", mse = 3)

data22 <- make_data(simulation37, simulation52, a_name = "Scenario 37 \n N = 5", b_name = "Scenario 52 \n N = 10",
                    xlab = "10 covariates, 10 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 3)

data23 <- make_data(simulation38, simulation53, a_name = "Scenario 38 \n N = 5", b_name = "Scenario 53 \n N = 10",
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 3)

data24 <- make_data(simulation39, simulation54, a_name = "Scenario 39 \n N = 5", b_name = "Scenario 54 \n N = 10",
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, small effect modification", mse = 3)

data25 <- make_data(simulation40, simulation55, a_name = "Scenario 40 \n N = 5", b_name = "Scenario 55 \n N = 10",
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, large effect modification", mse = 3)

data26 <- make_data(simulation41, simulation56, a_name = "Scenario 41 \n N = 5", b_name = "Scenario 56 \n N = 10",
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 3)

data27 <- make_data(simulation42, simulation57, a_name = "Scenario 42 \n N = 5", b_name = "Scenario 57 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 3)

data28 <- make_data(simulation43, simulation58, a_name = "Scenario 43 \n N = 5", b_name = "Scenario 58 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, small effect modification", mse = 3)

data29 <- make_data(simulation44, simulation59, a_name = "Scenario 44 \n N = 5", b_name = "Scenario 59 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, large effect modification", mse = 3)

data30 <- make_data(simulation45, simulation60, a_name = "Scenario 45 \n N = 5", b_name = "Scenario 60 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 3)




plot1 <- ggplot(data=data1$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data1$xlab)

plot2 <- ggplot(data=data2$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data2$xlab)

plot3 <- ggplot(data=data3$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data3$xlab)

plot4 <- ggplot(data=data4$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data4$xlab)

plot5 <- ggplot(data=data5$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data5$xlab)

plot6 <- ggplot(data=data6$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data6$xlab)

plot7 <- ggplot(data=data7$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data7$xlab)

plot8 <- ggplot(data=data8$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data8$xlab)

plot9 <- ggplot(data=data9$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data9$xlab)

plot10 <- ggplot(data=data10$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data10$xlab)

plot11 <- ggplot(data=data11$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data11$xlab)

plot12 <- ggplot(data=data12$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data12$xlab)

plot13 <- ggplot(data=data13$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data13$xlab)

plot13 <- ggplot(data=data13$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data13$xlab)

plot14 <- ggplot(data=data14$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data14$xlab)

plot15 <- ggplot(data=data15$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data15$xlab)

plot16 <- ggplot(data=data16$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data16$xlab)

plot17 <- ggplot(data=data17$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data17$xlab)

plot18 <- ggplot(data=data18$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data18$xlab)

plot19 <- ggplot(data=data19$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data19$xlab)

plot20 <- ggplot(data=data20$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data20$xlab)

plot21 <- ggplot(data=data21$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data21$xlab)

plot22 <- ggplot(data=data22$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data22$xlab)

plot23 <- ggplot(data=data23$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data23$xlab)

plot24 <- ggplot(data=data24$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data24$xlab)

plot25 <- ggplot(data=data25$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data25$xlab)

plot26 <- ggplot(data=data26$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data26$xlab)

plot27 <- ggplot(data=data27$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data27$xlab)

plot28 <- ggplot(data=data28$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data28$xlab)

plot29 <- ggplot(data=data29$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data29$xlab)

plot30 <- ggplot(data=data30$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data30$xlab)



grid.arrange(arrangeGrob(plot1,plot2, plot3,plot4,plot5,plot6,plot7,plot8,plot9, plot10, plot11, plot12, plot13, plot14, plot15,
                         top=fg, ncol=3, as.table = FALSE),
             arrangeGrob(plot16,plot17,plot18,plot19,plot20,plot21,plot22,plot23,plot24,plot25,plot26,plot27,plot28,plot29,plot30,
                         top=fg2, ncol=3, as.table = FALSE),
             left = textGrob("Treatment MSE", rot = 90, vjust = 0.5), ncol=2)





################ true effect modifiers


data1 <- make_data(simulation1, simulation16, a_name = "Scenario 1 \n N = 5", b_name = "Scenario 16 \n N = 10",
                   xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.2, no effect modification", mse = 2)

data2 <- make_data(simulation2, simulation17, a_name = "Scenario 2 \n N = 5", b_name = "Scenario 17 \n N = 10",
                   xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.5, no effect modification", mse = 2)

data3 <- make_data(simulation3, simulation18, a_name = "Scenario 3 \n N = 5", b_name = "Scenario 18 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, small effect modification", mse = 2)

data4 <- make_data(simulation4, simulation19, a_name = "Scenario 4 \n N = 5", b_name = "Scenario 19 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, small effect modification", mse = 2)

data5 <- make_data(simulation5, simulation20, a_name = "Scenario 5 \n N = 5", b_name = "Scenario 20 \n N = 10", 
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, large effect modification", mse = 2)

data6 <- make_data(simulation6, simulation21, a_name = "Scenario 6 \n N = 5", b_name = "Scenario 21 \n N = 10", 
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, large effect modification", mse = 2)

data7 <- make_data(simulation7, simulation22, a_name = "Scenario 7 \n N = 5", b_name = "Scenario 22 \n N = 10",
                   xlab = "10 covariates, 10 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 2)

data8 <- make_data(simulation8, simulation23, a_name = "Scenario 8 \n N = 5", b_name = "Scenario 23 \n N = 10", 
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 2)

data9 <- make_data(simulation9, simulation24, a_name = "Scenario 9 \n N = 5", b_name = "Scenario 24 \n N = 10", 
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, small effect modification", mse = 2)

data10 <- make_data(simulation10, simulation25, a_name = "Scenario 10 \n N = 5", b_name = "Scenario 25 \n N = 10", 
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, large effect modification", mse = 2)

data11 <- make_data(simulation11, simulation26, a_name = "Scenario 11 \n N = 5", b_name = "Scenario 26 \n N = 10", 
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 2)

data12 <- make_data(simulation12, simulation27, a_name = "Scenario 12 \n N = 5", b_name = "Scenario 27 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 2)

data13 <- make_data(simulation13, simulation28, a_name = "Scenario 13 \n N = 5", b_name = "Scenario 28 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, small effect modification", mse = 2)

data14 <- make_data(simulation14, simulation29, a_name = "Scenario 14 \n N = 5", b_name = "Scenario 29 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, large effect modification", mse = 2)

data15 <- make_data(simulation15, simulation30, a_name = "Scenario 15 \n N = 5", b_name = "Scenario 30 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 2)

## binary

data16 <- make_data(simulation31, simulation46, a_name = "Scenario 31 \n N = 5", b_name = "Scenario 46 \n N = 10",
                    xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.2, no effect modification", mse = 2)

data17 <- make_data(simulation32, simulation47, a_name = "Scenario 32 \n N = 5", b_name = "Scenario 47 \n N = 10",
                    xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.5, no effect modification", mse = 2)

data18 <- make_data(simulation33, simulation48, a_name = "Scenario 33 \n N = 5", b_name = "Scenario 48 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, small effect modification", mse = 2)

data19 <- make_data(simulation34, simulation49, a_name = "Scenario 34 \n N = 5", b_name = "Scenario 49 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, small effect modification", mse = 2)

data20 <- make_data(simulation35, simulation50, a_name = "Scenario 35 \n N = 5", b_name = "Scenario 50 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, large effect modification", mse = 2)

data21 <- make_data(simulation36, simulation51, a_name = "Scenario 36 \n N = 5", b_name = "Scenario 51 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, large effect modification", mse = 2)

data22 <- make_data(simulation37, simulation52, a_name = "Scenario 37 \n N = 5", b_name = "Scenario 52 \n N = 10",
                    xlab = "10 covariates, 10 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 2)

data23 <- make_data(simulation38, simulation53, a_name = "Scenario 38 \n N = 5", b_name = "Scenario 53 \n N = 10",
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 2)

data24 <- make_data(simulation39, simulation54, a_name = "Scenario 39 \n N = 5", b_name = "Scenario 54 \n N = 10",
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, small effect modification", mse = 2)

data25 <- make_data(simulation40, simulation55, a_name = "Scenario 40 \n N = 5", b_name = "Scenario 55 \n N = 10",
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, large effect modification", mse = 2)

data26 <- make_data(simulation41, simulation56, a_name = "Scenario 41 \n N = 5", b_name = "Scenario 56 \n N = 10",
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 2)

data27 <- make_data(simulation42, simulation57, a_name = "Scenario 42 \n N = 5", b_name = "Scenario 57 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 2)

data28 <- make_data(simulation43, simulation58, a_name = "Scenario 43 \n N = 5", b_name = "Scenario 58 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, small effect modification", mse = 2)

data29 <- make_data(simulation44, simulation59, a_name = "Scenario 44 \n N = 5", b_name = "Scenario 59 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, large effect modification", mse = 2)

data30 <- make_data(simulation45, simulation60, a_name = "Scenario 45 \n N = 5", b_name = "Scenario 60 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 2)




plot1 <- ggplot(data=data1$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data1$xlab)

plot2 <- ggplot(data=data2$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data2$xlab)

plot3 <- ggplot(data=data3$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data3$xlab)

plot4 <- ggplot(data=data4$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data4$xlab)

plot5 <- ggplot(data=data5$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data5$xlab)

plot6 <- ggplot(data=data6$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data6$xlab)

plot7 <- ggplot(data=data7$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data7$xlab)

plot8 <- ggplot(data=data8$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data8$xlab)

plot9 <- ggplot(data=data9$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data9$xlab)

plot10 <- ggplot(data=data10$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data10$xlab)

plot11 <- ggplot(data=data11$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data11$xlab)

plot12 <- ggplot(data=data12$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data12$xlab)

plot13 <- ggplot(data=data13$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data13$xlab)

plot13 <- ggplot(data=data13$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data13$xlab)

plot14 <- ggplot(data=data14$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data14$xlab)

plot15 <- ggplot(data=data15$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data15$xlab)

plot16 <- ggplot(data=data16$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data16$xlab)

plot17 <- ggplot(data=data17$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data17$xlab)

plot18 <- ggplot(data=data18$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data18$xlab)

plot19 <- ggplot(data=data19$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data19$xlab)

plot20 <- ggplot(data=data20$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data20$xlab)

plot21 <- ggplot(data=data21$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data21$xlab)

plot22 <- ggplot(data=data22$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data22$xlab)

plot23 <- ggplot(data=data23$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data23$xlab)

plot24 <- ggplot(data=data24$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data24$xlab)

plot25 <- ggplot(data=data25$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data25$xlab)

plot26 <- ggplot(data=data26$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data26$xlab)

plot27 <- ggplot(data=data27$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data27$xlab)

plot28 <- ggplot(data=data28$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data28$xlab)

plot29 <- ggplot(data=data29$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data29$xlab)

plot30 <- ggplot(data=data30$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data30$xlab)


grid.arrange(arrangeGrob(plot3,plot4,plot5,plot6,plot7,plot8,plot9, plot10, plot11, plot12, plot13, plot14, plot15,
                         top=fg, ncol=3, as.table = FALSE),
             arrangeGrob(plot18,plot19,plot20,plot21,plot22,plot23,plot24,plot25,plot26,plot27,plot28,plot29,plot30,
                         top=fg2, ncol=3, as.table = FALSE),
             left = textGrob("True effect modifier MSE", rot = 90, vjust = 0.5), ncol=2)





############## False effect modifiers


data1 <- make_data(simulation1, simulation16, a_name = "Scenario 1 \n N = 5", b_name = "Scenario 16 \n N = 10",
                   xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.2, no effect modification", mse = 1)

data2 <- make_data(simulation2, simulation17, a_name = "Scenario 2 \n N = 5", b_name = "Scenario 17 \n N = 10",
                   xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.5, no effect modification", mse = 1)

data3 <- make_data(simulation3, simulation18, a_name = "Scenario 3 \n N = 5", b_name = "Scenario 18 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, small effect modification", mse = 1)

data4 <- make_data(simulation4, simulation19, a_name = "Scenario 4 \n N = 5", b_name = "Scenario 19 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, small effect modification", mse = 1)

data5 <- make_data(simulation5, simulation20, a_name = "Scenario 5 \n N = 5", b_name = "Scenario 20 \n N = 10", 
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, large effect modification", mse = 1)

data6 <- make_data(simulation6, simulation21, a_name = "Scenario 6 \n N = 5", b_name = "Scenario 21 \n N = 10", 
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, large effect modification", mse = 1)

data7 <- make_data(simulation7, simulation22, a_name = "Scenario 7 \n N = 5", b_name = "Scenario 22 \n N = 10",
                   xlab = "10 covariates, 10 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 1)

data8 <- make_data(simulation8, simulation23, a_name = "Scenario 8 \n N = 5", b_name = "Scenario 23 \n N = 10", 
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 1)

data9 <- make_data(simulation9, simulation24, a_name = "Scenario 9 \n N = 5", b_name = "Scenario 24 \n N = 10", 
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, small effect modification", mse = 1)

data10 <- make_data(simulation10, simulation25, a_name = "Scenario 10 \n N = 5", b_name = "Scenario 25 \n N = 10", 
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, large effect modification", mse = 1)

data11 <- make_data(simulation11, simulation26, a_name = "Scenario 11 \n N = 5", b_name = "Scenario 26 \n N = 10", 
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 1)

data12 <- make_data(simulation12, simulation27, a_name = "Scenario 12 \n N = 5", b_name = "Scenario 27 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 1)

data13 <- make_data(simulation13, simulation28, a_name = "Scenario 13 \n N = 5", b_name = "Scenario 28 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, small effect modification", mse = 1)

data14 <- make_data(simulation14, simulation29, a_name = "Scenario 14 \n N = 5", b_name = "Scenario 29 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, large effect modification", mse = 1)

data15 <- make_data(simulation15, simulation30, a_name = "Scenario 15 \n N = 5", b_name = "Scenario 30 \n N = 10", 
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 1)

## binary

data16 <- make_data(simulation31, simulation46, a_name = "Scenario 31 \n N = 5", b_name = "Scenario 46 \n N = 10",
                    xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.2, no effect modification", mse = 1)

data17 <- make_data(simulation32, simulation47, a_name = "Scenario 32 \n N = 5", b_name = "Scenario 47 \n N = 10",
                    xlab = "10 covariates, 0 effect modifier \n \u03c4 = 0.5, no effect modification", mse = 1)

data18 <- make_data(simulation33, simulation48, a_name = "Scenario 33 \n N = 5", b_name = "Scenario 48 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, small effect modification", mse = 1)

data19 <- make_data(simulation34, simulation49, a_name = "Scenario 34 \n N = 5", b_name = "Scenario 49 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, small effect modification", mse = 1)

data20 <- make_data(simulation35, simulation50, a_name = "Scenario 35 \n N = 5", b_name = "Scenario 50 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.2, large effect modification", mse = 1)

data21 <- make_data(simulation36, simulation51, a_name = "Scenario 36 \n N = 5", b_name = "Scenario 51 \n N = 10",
                   xlab = "10 covariates, 1 effect modifier \n \u03c4 = 0.5, large effect modification", mse = 1)

data22 <- make_data(simulation37, simulation52, a_name = "Scenario 37 \n N = 5", b_name = "Scenario 52 \n N = 10",
                    xlab = "10 covariates, 10 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 1)

data23 <- make_data(simulation38, simulation53, a_name = "Scenario 38 \n N = 5", b_name = "Scenario 53 \n N = 10",
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 1)

data24 <- make_data(simulation39, simulation54, a_name = "Scenario 39 \n N = 5", b_name = "Scenario 54 \n N = 10",
                   xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, small effect modification", mse = 1)

data25 <- make_data(simulation40, simulation55, a_name = "Scenario 40 \n N = 5", b_name = "Scenario 55 \n N = 10",
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.2, large effect modification", mse = 1)

data26 <- make_data(simulation41, simulation56, a_name = "Scenario 41 \n N = 5", b_name = "Scenario 56 \n N = 10",
                    xlab = "15 covariates, 2 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 1)

data27 <- make_data(simulation42, simulation57, a_name = "Scenario 42 \n N = 5", b_name = "Scenario 57 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, small effect modification", mse = 1)

data28 <- make_data(simulation43, simulation58, a_name = "Scenario 43 \n N = 5", b_name = "Scenario 58 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, small effect modification", mse = 1)

data29 <- make_data(simulation44, simulation59, a_name = "Scenario 44 \n N = 5", b_name = "Scenario 59 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.2, large effect modification", mse = 1)

data30 <- make_data(simulation45, simulation60, a_name = "Scenario 45 \n N = 5", b_name = "Scenario 60 \n N = 10",
                    xlab = "15 covariates, 3 effect modifiers \n \u03c4 = 0.5, large effect modification", mse = 1)




plot1 <- ggplot(data=data1$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data1$xlab)

plot2 <- ggplot(data=data2$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data2$xlab)

plot3 <- ggplot(data=data3$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data3$xlab)

plot4 <- ggplot(data=data4$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data4$xlab)

plot5 <- ggplot(data=data5$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data5$xlab)

plot6 <- ggplot(data=data6$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data6$xlab)

plot7 <- ggplot(data=data7$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data7$xlab)

plot8 <- ggplot(data=data8$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data8$xlab)

plot9 <- ggplot(data=data9$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data9$xlab)

plot10 <- ggplot(data=data10$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data10$xlab)

plot11 <- ggplot(data=data11$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data11$xlab)

plot12 <- ggplot(data=data12$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data12$xlab)

plot13 <- ggplot(data=data13$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data13$xlab)

plot13 <- ggplot(data=data13$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data13$xlab)

plot14 <- ggplot(data=data14$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data14$xlab)

plot15 <- ggplot(data=data15$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data15$xlab)

plot16 <- ggplot(data=data16$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data16$xlab)

plot17 <- ggplot(data=data17$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data17$xlab)

plot18 <- ggplot(data=data18$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data18$xlab)

plot19 <- ggplot(data=data19$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data19$xlab)

plot20 <- ggplot(data=data20$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data20$xlab)

plot21 <- ggplot(data=data21$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data21$xlab)

plot22 <- ggplot(data=data22$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data22$xlab)

plot23 <- ggplot(data=data23$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data23$xlab)

plot24 <- ggplot(data=data24$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) +
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(data24$xlab)

plot25 <- ggplot(data=data25$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data25$xlab)

plot26 <- ggplot(data=data26$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data26$xlab)

plot27 <- ggplot(data=data27$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data27$xlab)

plot28 <- ggplot(data=data28$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data28$xlab)

plot29 <- ggplot(data=data29$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data29$xlab)

plot30 <- ggplot(data=data30$data, aes(x=models, y=error)) +
  geom_bar(stat="identity") + facet_grid(~simulations) + 
  theme(axis.title.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlab(data30$xlab)

grid.arrange(arrangeGrob(plot1,plot2,plot3,plot4,plot5,plot6,plot8,plot9, plot10, plot11, plot12, plot13, plot14, plot15,
                         top=fg, ncol=3, as.table = FALSE),
             arrangeGrob(plot16, plot17, plot18,plot19,plot20,plot21,plot23,plot24, plot25,plot26,plot27,plot28,plot29,plot30,
                         top=fg2, ncol=3, as.table = FALSE),
             left = textGrob("True effect modifier MSE", rot = 90, vjust = 0.5), ncol=2)
