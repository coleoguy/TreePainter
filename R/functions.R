treePaintR <- function(tree = NULL,
                       tip_states = NULL,
                       qmat = NULL,
                       iter = NULL,
                       rate.classes = NULL,
                       step = .5,
                       verbose = F){
  #TODO add some checks
  # is tree bifurcating
  # does qmat make sense
  # is iter sufficient

  # set rate classes
  steps <- floor(rate.classes/2)
  bot <- 1/(1 + steps * step)
  top <- 1 + steps * step
  rates <- c(seq(from = bot,
                 to = 1,
                 length.out = steps+1),
             seq(from = 1,
                 to = top,
                 length.out = steps+1)[-1])

  # add our new element to the tree
  tree$rates <- rep(ceiling(rate.classes/2),
                    Nedge(tree))
  # we will store the likelihood trace here
  lk.trace <- c()

  # get the starting likelihood
  lk1 <- asr_mk_model(tree = tree, tip_states = tip_states,
                      transition_matrix = qmat,
                      Nstates = ncol(qmat))$loglikelihood
  # this will be the primary loop
  for(j in 1:iter){
    if(verbose == T){
      if(j %% 500 == 0){
        cat(j, "\n")
      }
    }
    # this function tests a new rate only returns
    # a result if it is better
    result <- testNewRate(tree = tree, tip_states = tip_states,
                          transition_matrix = qmat,
                          Nstates = ncol(qmat),
                          rate.classes = rate.classes,
                          rates = rates, lk1 = lk1)
    if(length(result) > 0){
      tree <- result[[1]]
      lk1 <- result[[2]]
    }
    lk.trace[j] <- lk1
  }
  # here we prepare the results
  results <- list()
  class(tree) <- "rateTree"
  results[[1]] <- tree
  results[[2]] <- lk.trace
  results[[3]] <- rate.classes
  names(results) <- c("tree","lk.trace", "num.rates")
  return(results)
}

testNewRate <- function(tree = NULL, tip_states = NULL,
                        transition_matrix = NULL, Nstates = NULL,
                        rate.classes = NULL, rates = NULL, lk1 = NULL){
  # create a temp tree
  temp.tree <- tree
  # sample an edge at random
  edge <- sample(1:Nedge(tree), 1)
  # evaluate and draw the possible rates we could
  # assign to the sampled branch
  poss.rate <- getRates(tree = tree,
                       edge = edge,
                       rate.classes = rate.classes)
  # only bother doing something if it a change
  if(poss.rate != tree$rates[edge]){
    # apply the new rate
    temp.tree$rates[edge] <- poss.rate
    # change the branch lengths of all
    # branches based on their rates
    temp.tree$edge.length <- temp.tree$edge.length * rates[temp.tree$rates]
    # calculate the likelihood of the scaled tree
    lk2 <- asr_mk_model(tree = temp.tree,
                        tip_states = tip_states,
                        transition_matrix = transition_matrix,
                        Nstates = Nstates)$loglikelihood
    # see if the new rate is better
    if(lk2 > lk1){
      tree$rates <- temp.tree$rates
      lk1 <- lk2
      return(list(tree, lk1))
    }
  }else{
    return()
  }
}

# this function returns possible
# rates for a sampled branch
getRates <- function(tree = NULL, edge = NULL, rate.classes = NULL){
  # this gets the local rates
  local.rates <- getLocalRates(tree, edge, rate.classes)
  # this gets the intersection of possible rates and
  # returns the numeric that describes the rate class
  poss.rates <- getPossiRates(local.rates, rate.classes)
  return(poss.rates)
}

getLocalRates <- function(tree = NULL, edge = NULL, rate.classes){
  # get parent edge
  p.edge <- which(tree$edge[, 2] == tree$edge[edge, 1])
  # get daughter edge(s)
  d.edges <- which(tree$edge[, 1] == tree$edge[edge, 2])
  # if the branch has no parent we set it to the ML rate
  if(length(p.edge) == 0){
    rates <- median(1:rate.classes)
  }else{
    # get the rate of the parent edge
    rates <- tree$rates[p.edge]
  }
  # check to see if the branch has daughters
  if(length(d.edges) != 0){
    # collect the daugter rates
    rates <- c(rates, tree$rates[d.edges])
  }
  return(unique(rates))
}

getPossiRates <- function(local.rates = NULL, rate.classes = NULL){
  x <- c(-1, 0, 1)
  rate.mat <- as.data.frame(matrix(NA, 0, rate.classes))
  for(i in 1:length(local.rates)){
    z <- local.rates[i] + x
    z <- z[z > 0 & z <= rate.classes]
    rate.mat[i, z] <- 1
  }
  poss.rates <- as.numeric(which(colSums(rate.mat) == nrow(rate.mat)))
  if(length(poss.rates) > 1){
    poss.rates <- sample(poss.rates, 1)
  }
  return(poss.rates)
}

plotRateTree <- function(tree, rates,
                         scaled = F,
                         cols = NULL, bg = "white",
                         edge.width){
  foo <- par()
  if(scaled){
    tree$edge.length <- tree$edge.length * tree$rates
  }
  class(tree) <- "phylo"
  if(is.null(cols)){
    cols <- heat.colors(rates)[rates:1]
    cols[median(1:rates)] <- "darkgray"
  }
  if(bg != "white"){
    par(bg=bg)
  }
  plot(tree,
       show.tip.label = F,
       edge.color = cols[as.factor(tree$rates)],
       edge.width = edge.width)
  par(bg=foo$bg)
}
