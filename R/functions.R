treePaintR <- function(tree = NULL,
                       tip_states = NULL,
                       qmat = NULL,
                       iter = NULL,
                       rate.class = NULL,
                       rates = NULL,
                       verbose = T){
  lk.trace <- c()
  tree$rates <- rep(ceiling(length(rate.class)/2),
                    Nedge(tree))
  lk1 <- asr_mk_model(tree = tree, tip_states = tip_states,
                      transition_matrix = qmat,
                      Nstates = ncol(qmat))$loglikelihood
  for(j in 1:iter){
    if(j %% 500 == 0) cat(j, "\n")
    result <- testNewRate(tree = tree, tip_states = tip_states,
                          transition_matrix = qmat,
                          Nstates = ncol(qmat),
                          rate.class = rate.class,
                          rates = rates, lk1 = lk1)
    if(length(result) > 0){
      tree <- result[[1]]
      lk1 <- result[[2]]
    }
    lk.trace[j] <- lk1
  }
  results <- list()
  results[[1]] <- tree
  results[[2]] <- lk.trace
  names(results) <- c("tree","lk.trace")
  return(results)
}
testNewRate <- function(tree = NULL, tip_states = NULL,
                        transition_matrix = NULL, Nstates = NULL,
                        rate.class = NULL, rates = NULL, lk1 = NULL){
  temp.tree <- tree
  edge <- sample(1:Nedge(tree), 1)
  poss.rate <- getRates(tree = tree,
                       edge = edge,
                       rate.class = rate.class)
  if(poss.rate != tree$rates[edge]){
    temp.tree$rates[edge] <- poss.rate
    temp.tree$edge.length <- temp.tree$edge.length * rates[temp.tree$rates]
    lk2 <- asr_mk_model(tree = temp.tree,
                        tip_states = tip_states,
                        transition_matrix = transition_matrix,
                        Nstates = Nstates)$loglikelihood
    if(lk2 > lk1){
      tree$rates <- temp.tree$rates
      lk1 <- lk2
      return(list(tree, lk1))
    }
  }else{
    return()
  }
}

getRates <- function(tree = NULL, edge = NULL, rate.class = NULL){
  local.rates <- getLocalRates(tree, edge)
  poss.rates <- getPossiRates(local.rates, rate.class)
  return(poss.rates)
}
getLocalRates <- function(tree = NULL, edge = NULL){
  p.edge <- which(tree$edge[, 2] == tree$edge[edge, 1])
  d.edges <- which(tree$edge[, 1] == tree$edge[edge, 2])
  if(length(p.edge) == 0){
    rates <- median(rate.class)
  }else{
    rates <- tree$rates[p.edge]
  }
  if(length(d.edges) != 0){
    rates <- c(rates, tree$rates[d.edges])
  }
  return(unique(rates))
}
getPossiRates <- function(local.rates = NULL, rate.class = NULL){
  x <- c(-1, 0, 1)
  rate.mat <- as.data.frame(matrix(NA, 0, length(rate.class)))
  for(i in 1:length(local.rates)){
    z <- local.rates[i] + x
    z <- z[z > 0 & z <= length(rate.class)]
    rate.mat[i, z] <- 1
  }
  poss.rates <- as.numeric(which(colSums(rate.mat) == nrow(rate.mat)))
  if(length(poss.rates) > 1){
    poss.rates <- sample(poss.rates, 1)
  }
  return(poss.rates)
}


