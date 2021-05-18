library(diversitree)
library(castor)
source("functions.R")

qmat <- matrix(c(-.3,.3,.3,-.3), 2,2)
nTrees <- 1
trees <- traits <- list()
fast.branches <- list()
ntips <- c(200)
for(i in 1:nTrees){
  trees[[i]] <- trees(pars=c(3,1),
                      type="bd",
                      max.taxa = ntips[i])[[1]]
  # scale tree to unit length
 # depth <- max(branching.times(trees[[i]]))
#  trees[[i]]$edge.length <- trees[[i]]$edge.length / depth
  working <- T
  counter <- 0
  while(working){
    print(counter + 1)
    counter <- counter + 1
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

sim.tree <- trees[[i]]
sim.tree$edge.length[fast.branches[[i]]] <- sim.tree$edge.length[fast.branches[[i]]] * 5
# traits[[i]] <- as.vector(sim.character(tree = sim.tree,
#                                        pars = qmat,
#                                        model = "mkn",x0 = 2))

traits[[i]] <- simulate_mk_model(tree = sim.tree, Q = qmat, root_probabilities = c(0,1),include_nodes = F,Nsimulations = 1,drop_dims = T)$tip_states

# simulate_mk_model(tree = sim.tree, Q = qmat, root_probabilities = c(0,1),include_nodes = F,Nsimulations = 1,drop_dims = T)$tip_states

# make sure that two traits are present
if(length(unique(traits[[i]])) == 1){
  working <- T
  while (working) {
    # traits[[i]] <- as.vector(sim.character(tree = sim.tree,
    #                                        pars = qmat,
    #                                        model = "mkn",x0 = 2))
    
    traits[[i]] <- simulate_mk_model(tree = sim.tree, Q = qmat, root_probabilities = c(0,1),include_nodes = F,Nsimulations = 1,drop_dims = T)$tip_states
    if(length(unique(traits[[i]])) == 2){
      working <- F
    }
  }
}

plot(sim.tree, show.tip.label = F)
tiplabels(pch = 16, col = c("red", "blue")[traits[[1]]], cex = .5, offset = .01)

plot(x = NULL, y = NULL,
     xlim = c(1,10000),
     ylim = c(-90,-70))


for(j in 1:20){
x <- treePaintR(tree = trees[[i]],
           tip_states = traits[[i]],
           qmat = qmat,
           iter = 500,
           rate.classes = 5,
           step = 0.1,iter.check = T,
           iter.check.interval = 1000)
lines(x$lk.trace, col = rainbow(20)[j])
}




plot.rateTree(x$tree, rates = x$num.rates,edge.width = 1, scaled = F)
hist(x$tree$rates)

# true positives
th <- sum(x$tree$rates[fast.branches[[i]]]> median(1:x$num.rates)) / length(fast.branches[[i]])
tl <- sum(x$tree$rates[-fast.branches[[i]]] < median(1:x$num.rates)) / length(x$tree$rates[-fast.branches[[i]]])
# false positives
fh <- sum(x$tree$rates[-fast.branches[[i]]] > median(1:x$num.rates)) / length(x$tree$rates[-fast.branches[[i]]])
fl <- sum(x$tree$rates[fast.branches[[i]]] < median(1:x$num.rates)) / length(fast.branches[[i]])

th
tl
fh
fl


col <- rep("black",Nedge(trees[[1]]))
col[x$tree$rates > median(1:x$num.rates)] <- "red"
col[x$tree$rates == median(1:x$num.rates)] <- "black"
col[x$tree$rates < median(1:x$num.rates)] <- "blue"

plot(trees[[1]], edge.color = col, show.tip.label = F)
tiplabels(pch = 16, col = c("red", "blue")[traits[[1]]], cex = .5, offset = .01)
plot(sim.tree, show.tip.label = F)
tiplabels(pch = 16, col = c("red", "blue")[traits[[1]]], cex = .5, offset = .01)
RC <- 7
SS <- 0.5

{
  steps <- floor(RC/2)
bot <- 1/(1 + steps * SS)
top <- 1 + steps * SS
rates <- c(seq(from = bot,
               to = 1,
               length.out = steps+1),
           seq(from = 1,
               to = top,
               length.out = steps+1)[-1])
}
rates
plot(rates)
