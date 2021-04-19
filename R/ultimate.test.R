library(diversitree)
source("functions.R")
qmat <- matrix(c(-.3,.3,.3,-.3), 2,2)
trees <- traits <- list()
fast.branches <- list()
ntips <- c(50, 100, 150, 200)
for(i in 1:4){
  trees[[i]] <- trees(pars=c(3,1),
                      type="bd",
                      max.taxa = ntips[i])[[1]]
  # scale tree to unit length
  depth <- max(branching.times(trees[[i]]))
  trees[[i]]$edge.length <- trees[[i]]$edge.length / depth
  working <- T
  while(working){
    hit <- sample(trees[[i]]$edge[,1], 1)
    tips <- get.descendants(node = hit,
                            tree = trees[[i]], tips.only = T)
    branches <- get.descendants(node = hit,
                                tree = trees[[i]], edge.index = T)
    if(length(tips) > .2*length(trees[[i]]$tip.label) &
       length(tips) < .8*length(trees[[i]]$tip.label)){
      working <- F
      fast.branches[[i]] <- branches
    }
  }
}
library(castor)

# test in model performance / over fitting
rateClasses <- seq(from = 3, to = 21, by = 2)
stepSize <- seq(from = 1, to = 30, by = 2)
reps <- 5
scaler <- c(2,5,25,100)

# # get a single tree and traits
# tree <- trees[[i]]
# trait <- traits[[i]]

dat <- as.data.frame(matrix(data = NA,
                            nrow = length(trees) * reps * length(rateClasses) * length(stepSize) * length(scaler),
                            ncol = 8))

colnames(dat) <- c("tree", "rep", "Ntips", "rateClasses", "stepSize", "scaler", "truePositives", "falsePositives")

j <- 1
for(m in 1:length(trees)){
  for(i in 1:reps){
    for(k in 1:length(rateClasses)){
      for(l in 1:length(stepSize)){
        for(n in 1: length(scaler)){
          sim.tree <- trees[[m]]
          sim.tree$edge.length[branches] <- sim.tree$edge.length[branches] * scaler[n]
          traits[[m]] <- as.vector(sim.character(tree = sim.tree,
                                                 pars = qmat,
                                                 model = "mkn",x0 = 2))
          if(length(unique(traits[[m]])) == 1){
            working <- T
            while (working) {
              traits[[m]] <- as.vector(sim.character(tree = sim.tree,
                                                     pars = qmat,
                                                     model = "mkn",x0 = 2))
              if(length(unique(traits[[m]])) == 2){
                working <- F
              }
            }
          }
          #fill the data table
          dat$tree[j] <- paste("tree", m)
          dat$rep[j] <- paste("rep", i)
          dat$Ntips[j] <- Ntip(trees[[m]])
          dat$rateClasses[j] <- rateClasses[k]
          dat$stepSize[j] <- stepSize[l]
          dat$scaler[j] <- scaler[n]
          # print iteration
          print(paste("iteration", j))
          # fit model
          fit <- treePaintR(tree = trees[[m]],
                            tip_states = traits[[m]],
                            qmat = qmat,
                            iter = 6000,
                            rate.classes = rateClasses[k],
                            step = stepSize[l])
          # true positives
          tp <- sum(fit$tree$rates[fast.branches[[m]]] > median(1:fit$num.rates)) / length(fast.branches[[m]])
          # false positives
          fp <- sum(fit$tree$rates[-fast.branches[[m]]] > median(1:fit$num.rates)) / length(fit$tree$rates[-fast.branches[[m]]])
          dat$truePositives[j] <- tp
          dat$falsePositives[j] <- fp
          # counter
          j <- j+1
        }
      }
    }
  }
  print(paste("rep", i, "complete"))
}

par(mfcol = c(2,2))

# plot
## rate classes vs true positives
plot(x = jitter(dat$rateClasses[dat$tree == paste("tree", i)],.6),
     y = jitter(dat$truePositivesdat$tree == paste("tree", i),2),
     xlab = "Rate class",
     ylab = "True positive",
     pch = 16,
     col = rainbow(reps,alpha = .5)[as.factor(dat$rep[dat$tree == paste("tree", i)])])

## rate classes vs false positives
plot(x = jitter(dat$rateClasses[dat$tree == paste("tree", i)], .6),
     y = jitter(dat$falsePositives[dat$tree == paste("tree", i)],2),
     xlab = "Rate class",
     ylab = "False positive",
     pch = 16,
     col = rainbow(reps,alpha = .5)[as.factor(dat$rep[dat$tree == paste("tree", i)])])

## step size vs true positives
plot(x = jitter(dat$stepSize[dat$tree == paste("tree", i)], .6),
     y = jitter(dat$truePositives[dat$tree == paste("tree", i)], 2),
     xlab = "Step size",
     ylab = "True positive",
     pch = 16,
     col = rainbow(reps,alpha = .5)[as.factor(dat$rep[dat$tree == paste("tree", i)])])

## step size vs false positives
plot(x = jitter(dat$stepSize[dat$tree == paste("tree", i)], .6),
     y = jitter(dat$falsePositives[dat$tree == paste("tree", i)], 2),
     xlab = "Step size",
     ylab = "False positive",
     pch = 16,
     col = rainbow(reps,alpha = .5)[as.factor(dat$rep[dat$tree == paste("tree", i)])],
     cex = 1)

fit <- treePaintR(tree = tree,
                  tip_states = trait,
                  qmat = qmat,
                  iter = 10000,
                  rate.classes = dat$rateClasses[which.max(dat$truePositives)],
                  step = dat$stepSize[which.max(dat$truePositives)])

par(mfcol = c(1,2))
plot(tree = fit[[1]],
     rates = fit$num.rates,
     scaled = F,
     cols=NULL,
     bg="lightgray", edge.width = 1)
plot(sim.tree, show.tip.label = F)
tiplabels(pch = 16, col = c("red", "blue")[traits[[1]]], cex = .5, offset = .01)




for(i in 1:4){
  # plot
  ## rate classes vs true positives
  plot(x = jitter(dat$rateClasses[dat$tree == paste("tree", i)],.6),
       y = jitter(dat$truePositives[dat$tree == paste("tree", i)],2),
       xlab = "Rate class",
       ylab = "True positive",
       pch = 16,
       col = rainbow(reps,alpha = .5)[as.factor(dat$rep[dat$tree == paste("tree", i)])])
  
  ## rate classes vs false positives
  plot(x = jitter(dat$rateClasses[dat$tree == paste("tree", i)], .6),
       y = jitter(dat$falsePositives[dat$tree == paste("tree", i)],2),
       xlab = "Rate class",
       ylab = "False positive",
       pch = 16,
       col = rainbow(reps,alpha = .5)[as.factor(dat$rep[dat$tree == paste("tree", i)])])
  
  ## step size vs true positives
  plot(x = jitter(dat$stepSize[dat$tree == paste("tree", i)], .6),
       y = jitter(dat$truePositives[dat$tree == paste("tree", i)], 2),
       xlab = "Step size",
       ylab = "True positive",
       pch = 16,
       col = rainbow(reps,alpha = .5)[as.factor(dat$rep[dat$tree == paste("tree", i)])])
  
  ## step size vs false positives
  plot(x = jitter(dat$stepSize[dat$tree == paste("tree", i)], .6),
       y = jitter(dat$falsePositives[dat$tree == paste("tree", i)], 2),
       xlab = "Step size",
       ylab = "False positive",
       pch = 16,
       col = rainbow(reps,alpha = .5)[as.factor(dat$rep[dat$tree == paste("tree", i)])],
       cex = 1)
}

library(ggplot2)
#tp vs rc
ggplot(dat, aes(y=truePositives, x=rateClasses)) + geom_point(aes(colour=as.factor(Ntips)), stat="identity", position="jitter", alpha=0.7, size=1) + facet_grid(rep ~ scaler) + theme_grey() + theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + scale_size(range=c(1, 2)) + guides(colour=guide_legend(title="Ntips")) + xlab("rateClasses") + ylab("truePositives")

#fp vs rc
ggplot(dat, aes(y=falsePositives, x=rateClasses)) + geom_point(aes(colour=as.factor(Ntips)), stat="identity", position="jitter", alpha=0.7, size=1) + facet_grid(rep ~ scaler) + theme_grey() + theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + scale_size(range=c(1, 2)) + guides(colour=guide_legend(title="Ntips")) + xlab("rateClasses") + ylab("falsePositives")

# tp vs ss
ggplot(dat, aes(y=truePositives, x=stepSize)) + geom_point(aes(colour=as.factor(Ntips)), stat="identity", position="jitter", alpha=0.7, size=1) + facet_grid(rep ~ scaler) + theme_grey() + theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + scale_size(range=c(1, 2)) + guides(colour=guide_legend(title="Ntips")) + xlab("stepSize") + ylab("truePositives")

# fp vs ss
ggplot(dat, aes(y=falsePositives, x=stepSize)) + geom_point(aes(colour=as.factor(Ntips)), stat="identity", position="jitter", alpha=0.7, size=1) + facet_grid(rep ~ scaler) + theme_grey() + theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + scale_size(range=c(1, 2)) + guides(colour=guide_legend(title="Ntips")) + xlab("stepSize") + ylab("falsePositives")

     