library(diversitree)
source("functions.R")
qmat <- matrix(c(-.3,.3,.3,-.3), 2,2)
trees <- traits <- list()
fast.branches <- list()
for(i in 101){
  trees[[i]] <- trees(pars=c(3,1),
                      type="bd",
                      max.taxa = (i+99))[[1]]
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
  sim.tree <- trees[[i]]
  sim.tree$edge.length[branches] <- sim.tree$edge.length[branches] * 5
  traits[[i]] <- as.vector(sim.character(tree = sim.tree,
                                         pars = qmat,
                                         model = "mkn",x0 = 2))
}
library(castor)

# test in model performance / over fitting
rateClasses <- seq(from = 3, to = 21, by = 2)
stepSize <- seq(from = 1, to = 30, by = 2)
reps <- 5
# scaler <- 2:5

# get a single tree and traits
tree <- trees[[i]]
trait <- traits[[i]]

dat <- as.data.frame(matrix(data = NA,
                            nrow = reps * length(rateClasses) * length(stepSize),
                            ncol = 6))

colnames(dat) <- c("rep", "Ntips", "rateClasses", "stepSize", "truePositives", "falsePositives")

j <- 1

for(i in 1:reps){
  for(k in 1:length(rateClasses)){
    for(l in 1:length(stepSize)){
      #fill the data table
      dat$rep[j] <- paste("rep", i)
      dat$Ntips[j] <- Ntip(tree)
      dat$rateClasses[j] <- rateClasses[k]
      dat$stepSize[j] <- stepSize[l]
      # print iteration
      print(paste("iteration", j))
      # fit model
      fit <- treePaintR(tree = tree,
                        tip_states = trait,
                        qmat = qmat,
                        iter = 10000,
                        rate.classes = rateClasses[k],
                        step = stepSize[l])
      # true positives
      tp <- sum(fit$tree$rates[fast.branches[[101]]] > median(1:fit$num.rates)) / length(fast.branches[[101]])
      # false positives
      fp <- sum(fit$tree$rates[-fast.branches[[101]]] > median(1:fit$num.rates)) / length(fit$tree$rates[-fast.branches[[101]]])
      dat$truePositives[j] <- tp
      dat$falsePositives[j] <- fp
      # counter
      j <- j+1
    }
  }
  print(paste("rep", i, "complete"))
}

par(mfcol = c(2,2))

# plot
## rate classes vs true positives
plot(x = jitter(dat$rateClasses,.6),
     y = jitter(dat$truePositives,2),
     xlab = "Rate class",
     ylab = "True positive",
     pch = 16,
     col = rainbow(reps,alpha = .5)[as.factor(dat$rep)])

## rate classes vs false positives
plot(x = jitter(dat$rateClasses, .6),
     y = jitter(dat$falsePositives,2),
     xlab = "Rate class",
     ylab = "False positive",
     pch = 16,
     col = rainbow(reps,alpha = .5)[as.factor(dat$rep)])

## step size vs true positives
plot(x = jitter(dat$stepSize, .6),
     y = jitter(dat$truePositives, 2),
     xlab = "Step size",
     ylab = "True positive",
     pch = 16,
     col = rainbow(reps,alpha = .5)[as.factor(dat$rep)])

## step size vs false positives
plot(x = jitter(dat$stepSize, .6),
     y = jitter(dat$falsePositives, 2),
     xlab = "Step size",
     ylab = "False positive",
     pch = 16,
     col = rainbow(reps,alpha = .5)[as.factor(dat$rep)],
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
tiplabels(pch = 16, col = c("red", "blue")[traits[[101]]], cex = .5, offset = .01)


