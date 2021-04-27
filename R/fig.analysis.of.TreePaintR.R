# load libraries
library(plotly)
#load data 
load("ultimate.test.RData")
# remove unwanted data
rm(list = ls()[-3])
# set pars for subsetting
Ntips <- 200
scaler <- 5
# scaler
StepSize <- dat$stepSize[dat$Ntips == Ntips & dat$rep == "rep 1" & dat$scaler == scaler]
# rate classes
RateClass <- dat$rateClasses[dat$Ntips == Ntips & dat$rep == "rep 1" & dat$scaler == scaler]
# true positives
TruePositives <- rowMeans(data.frame(list(dat$truePositives[dat$Ntips == Ntips & dat$rep == "rep 1" & dat$scaler == scaler],
                                 dat$truePositives[dat$Ntips == Ntips & dat$rep == "rep 2" & dat$scaler == scaler],
                                 dat$truePositives[dat$Ntips == Ntips & dat$rep == "rep 3" & dat$scaler == scaler],
                                 dat$truePositives[dat$Ntips == Ntips & dat$rep == "rep 4" & dat$scaler == scaler],
                                 dat$truePositives[dat$Ntips == Ntips & dat$rep == "rep 5" & dat$scaler == scaler])))
# false positives
FalsePositives <- rowMeans(data.frame(list(dat$falsePositives[dat$Ntips == Ntips & dat$rep == "rep 1" & dat$scaler == scaler],
                                 dat$falsePositives[dat$Ntips == Ntips & dat$rep == "rep 2" & dat$scaler == scaler],
                                 dat$falsePositives[dat$Ntips == Ntips & dat$rep == "rep 3" & dat$scaler == scaler],
                                 dat$falsePositives[dat$Ntips == Ntips & dat$rep == "rep 4" & dat$scaler == scaler],
                                 dat$falsePositives[dat$Ntips == Ntips & dat$rep == "rep 5" & dat$scaler == scaler])))
# get figures
fig1 <- plot_ly(x = ~RateClass,
        y = ~StepSize,
        z = ~TruePositives, 
        type = "contour")
fig2 <- plot_ly(x = ~RateClass,
        y = ~StepSize,
        z = ~FalsePositives, 
        type = "contour")
# combine figures
fig <- subplot(fig1,fig2, shareX = T,shareY = T, titleX = T)
# plot
fig










