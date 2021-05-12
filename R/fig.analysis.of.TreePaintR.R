
#load data 
load("analysis.of.TreePaintR.scaler.5.RData")
# remove unwanted data
rm(list = ls()[-3])


StepSize <- unique(dat$stepSize)
RateClass <- unique(dat$rateClasses)

TruePositives <- FalsePositives <- matrix(data = NA,
                                          ncol = length(StepSize),
                                          nrow = length(RateClass))

colnames(TruePositives) <- colnames(FalsePositives) <- StepSize
rownames(TruePositives) <- rownames(FalsePositives) <- RateClass

for(i in 1:length(RateClass)){
        for(j in 1:length(StepSize)){
                TruePositives[i,j] <- mean(dat$truePositives[dat$rateClasses == RateClass[i] & dat$stepSize == StepSize[j]])
                FalsePositives[i,j] <- mean(dat$falsePositives[dat$rateClasses == RateClass[i] & dat$stepSize == StepSize[j]])
        }
}

filled.contour(x = RateClass, 
               y = StepSize, 
               z = TruePositives * 100,
               plot.title = title(main = "True positives",
                                  xlab = "Rate Class", ylab = "Step size"),
               key.title = title(main = "% of\n True positives",cex.main = .7),
               color = function(n) hcl.colors(n, "viridis"),
               levels = pretty(range(TruePositives * 100, FalsePositives * 100),20))

filled.contour(x = RateClass, 
               y = StepSize, 
               z = FalsePositives * 100, 
               plot.title = title(main = "False positives",
                                  xlab = "Rate Class", ylab = "Step size"),
               key.title = title(main = "% of\n False positives",cex.main = .7),
               color = function(n) hcl.colors(n, "viridis"),
               levels = pretty(range(TruePositives * 100, FalsePositives * 100),20))




# load libraries
library(plotly)
# set pars for subsetting
Ntips <- 200
scaler <- 5
# scaler
StepSize <- dat$stepSize[dat$Ntips == Ntips & dat$rep == "rep 1" & dat$scaler == scaler]
# rate classes
RateClass <- dat$rateClasses[dat$Ntips == Ntips & dat$rep == "rep 1" & dat$scaler == scaler]
# true positives and false positives
x <- y <- list()
for(i in 1:50){
        x[[i]] <- dat$truePositives[dat$Ntips == Ntips & dat$rep == paste("rep", i) & dat$scaler == scaler]
        y[[i]] <- dat$falsePositives[dat$Ntips == Ntips & dat$rep == paste("rep", i) & dat$scaler == scaler]
}
TruePositives <- rowMeans(data.frame(x))
FalsePositives <- rowMeans(data.frame(y))

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