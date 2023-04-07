
# Given any particular model, construct the proposal model for
# the grouped Metropolis Hastings algorithm. At any given stage,
# We can either add, remove or exchange any other variables, the 
# selection will be done proportion to weights specified by the 
# initial multivariate group lasso estimates

proposal <- function(addvar, cur.model, n_covar, weights){
  if (addvar == "Add" | length(cur.model) == 0){
    VarNotUsed   <- setdiff(1:n_covar, cur.model) # Variables not in the current model
    # if all the current model is the full model
    if (length(VarNotUsed) == 0){
      new.model  <- cur.model
      qratio     <- 1
    }else{
      prob       <- weights[VarNotUsed]/sum(weights[VarNotUsed])
      if (length(prob)==1){
        newvar     <- VarNotUsed
      }else{
        newvar     <- sample(VarNotUsed,1,prob = prob)
      }
      new.model  <- sort(c(cur.model, newvar))
      qratio     <- ( ((1 / weights[newvar]) / sum(1 / weights[new.model])) / (weights[newvar] / sum(weights[VarNotUsed])) )
    }
  }
  else if (addvar == "Remove"){
    if (length(cur.model) == 1){
      new.model  <- cur.model
      qratio     <- 1
    }else{
      invprob    <- ( 1 / weights[cur.model] ) / sum( 1 / weights[cur.model] )
      discardvar <- sample(cur.model,1,prob = invprob)
      new.model  <- setdiff(cur.model, discardvar)
      newVarNotUsed <- setdiff(1:n_covar, new.model)
      qratio     <- ( (weights[discardvar] / sum(weights[newVarNotUsed])) / ((1 / weights[discardvar]) / sum(1 / weights[cur.model])) )
    }
  }
  else{
    VarNotUsed   <- setdiff(1:n_covar, cur.model) # Variables not in the current model
    # if all the current model is the full model
    if (length(VarNotUsed) == 0){
      new.model  <- cur.model
      qratio     <- 1
    }else{
      prob       <- weights[VarNotUsed]/sum(weights[VarNotUsed])
      if (length(prob)==1){
        newvar   <- VarNotUsed
      }else{
        newvar     <- sample(VarNotUsed,1,prob = prob)
      }
      invprob    <- ( 1 / weights[cur.model] ) / sum( 1 / weights[cur.model] )
      if (length(invprob)==1){
        discardvar <- cur.model
      }else{
        discardvar <- sample(cur.model,1,prob = invprob)
      }
      new.model  <- sort(setdiff(c(cur.model, newvar),discardvar))
      newVarNotUsed <- setdiff(1:n_covar, new.model)
      qratio     <-  (weights[discardvar] / sum(weights[newVarNotUsed]))*((1 / weights[newvar]) / sum(1 / weights[new.model])) / (((1 / weights[discardvar]) / sum(1 / weights[cur.model]))*(weights[newvar] / sum(weights[VarNotUsed]))) 
    }
  }
  list("newmodel"=new.model, "qratio"=qratio)
}