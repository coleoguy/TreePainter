treePaintR <- function(tree = NULL,
                       tip_states = NULL,
                       qmat = NULL,
                       iter = NULL,
                       rate.classes = NULL,
                       step = .5,
                       verbose = F,
                       iter.check = NULL,
                       iter.check.interval = NULL){
  # do the checks
  # add checks
  # checks related to the tree
  # presence of tree
  if(is.null(tree)){
    stop("No imput tree")
  }else{
    # class of tree
    if(class(tree) != "phylo"){
      stop("tree should class of 'Phylo'")
    }
    # bifurcation of tree
    if(!is.binary(tree)){
      stop("tree should be stricly bifurcating")
    }
    # rootedness
    if(!is.rooted(tree)){
      stop("tree should be rooted")
    }
  }
  # checks related to tip states
  # presence of tip states
  if(is.null(tip_states)){
    stop("tip states are missing")
  }else{
    # length of tip states should be equal to the number of tips
    if(length(tip_states) != Ntip(tree)){
      stop("length of tip states does not match with the number of tips")
    }
    # tip states should be integer
    if(!is.numeric(tip_states)){
      stop("tip states should integers")
    }
    # tip states should start from 1
    if(min(range(tip_states)) != 1){
      stop("Starting value of tip states should be 1")
    }
  }
  # checks related to qmat
  # presence of qmat
  if(is.null(qmat)){
    stop("missing Q matrix")
  }else{
    if(nrow(qmat) != ncol(qmat)){
      stop("number of rows and columns in the Q Matrix does not match")
    }
    if(any(rowSums(qmat) != 0)){
      stop(paste("row(s)", which(rowSums(qmat) != 0), "does not add up to zero"))
    }
  }
  # checks related to iterations
  # presence of iterations
  if(is.null(iter)){
    stop("iter is missing")
  }else{
    # iter should be single value
    if(length(iter) > 1){
      stop("iter should be single numeric value")
    }
    # iter should be positive
    if(iter < 0){
      stop("iter should be positive")
    }
    # iter should be numeric
    if(!is.numeric(iter)){
      stop("iter should be of class an integer")
    }
    # iter should not contain decimal values
    if(iter %% 1 != 0){
      stop("rate classes should not contain decimal values")
    }
  }
  # checks related to rate class
  # presence of rate class
  if(is.null(rate.classes)){
    stop("rate classes are missing")
  }else{
    # rate class should be single value
    if(length(rate.classes) > 1){
      stop("rate classes should be single numeric value")
    }
    # rate class should be positive
    if(rate.classes < 0){
      stop("rate class should be positive")
    }
    # rate class should be numeric
    if(!is.numeric(rate.classes)){
      stop("rate classes should be an integer")
    }
    # rate class should not contain decimal values
    if(rate.classes %% 1 != 0){
      stop("rate classes should not contain decimal values")
    }
  }
  # checks related to iter check
  if(is.null(iter.check)){
    iter.check <- F
  }else{
    if(!(is.logical(iter.check))){
      stop("iter.check should be logical")
    }
  }
  # checks related to iter check interval
  if(is.null(iter.check.interval)){
    iter.check.interval <- 1000
  }else{
    if(iter.check.interval %% 1 != 0){
      print("iter.check.interval is not an integer. Setting it to the default value of 1000 generations")
      iter.check.interval <- 1000
    }
  }
  # iter and iter check interval cannot be same values
  if(iter.check == T){
    if(iter == iter.check.interval){
      stop("iter and iter check interval has to be different values. However the difference can be as small as 1")
    }
  }
  
  # set rate classes
  #TODO make this step better
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
  
  if(iter.check == T){
    if(iter.check.interval > iter){
      deviation <- sd(lk.trace[(iter - (iter * 0.5)):iter])      
    }else{
      deviation <- sd(lk.trace[(iter - iter.check.interval):iter])
    }
    if(deviation != 0){
      print(paste("Performing additional iterations untill no increase in the likelihood is seen within", iter.check.interval, "is seen."))
      counter <- T
      while (counter) {
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
        j <- j + 1
        if(iter.check.interval > length(lk.trace)){
          deviation <- sd(lk.trace[(length(lk.trace) - (length(lk.trace) *  (length(lk.trace)/iter.check.interval))):length(lk.trace)])
        }else{
          deviation <- sd(lk.trace[(length(lk.trace) - iter.check.interval):length(lk.trace)])
        }
        if(deviation == 0){
          counter <- F
        }
      }
    }
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

# this function returns possible
# rates of local branches
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

# this function returns possible
# rates for a set of local rates
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

# this function plots the tree painted/scaled by rate
plot.rateTree <- function(tree, rates,
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
