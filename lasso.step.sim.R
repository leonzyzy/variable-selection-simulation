library(mc2d)
library(MASS)
library(mvtnorm)
library(Matrix)
library(magic)
library(mixmeta)
library(rlist)

#====#====#====#====#====#==== function para define.
beta = matrix(rep(c(3.8, -2.8, 3.5, 0, -0.8, 2.2, -3.2, -4.2, 0, 1.8, 0.3, 1.1, -0.6, 0),
                  c(1, 1, 1, 4, 1, 1, 1, 1, 10, 1, 1, 1, 1, 5)), nrow = 30)

block.group = c(3,4,4,3,4,3,4,3,2) 
block.num = 1:30
nrow = 30
ncol = 30
r = c(0.25, 0.50, 0.75)


#====#====#====#====#====#========================= sim data function
#### This function allows you to input different correlation r, block matrix size and group.
sim_data = function(r, vec, group, nrow, ncol, beta, seed){
  ##define cov matrix
  cov.matrix = matrix(nrow=nrow, ncol=ncol)
  
  #fit values using r^|i-j|
  for(i in 1:30){
    for(j in 1:30){
      cov.matrix[i,j] = r^(abs(i-j))
    }
  }
  #create a list
  block.matrix.list = list()
  div.list = list()
  div.list[[1]] = vec[1:group[1]]
  
  #use loop to div number into values sum. i.e 30 = 3+4+4+3+4+3+4+3+2
  for(i in 2:length(group)){
    left.indx = div.list[[i-1]][length(div.list[[i-1]])] + 1
    right.indx = sum(group[1:i])
    
    div.list[[i]] = vec[(left.indx):(right.indx)]
  }
  
  #use loop create a block matrix based on different size.
  for(i in 1:length(group)){
    block.matrix.list[[i]] = cov.matrix[div.list[[i]], div.list[[i]]]
  }
  
  #create a block matrix
  block.matrix = bdiagMat(block.matrix.list)
  
  #simulate X with 100X30 dim from muti-normal.
  set.seed(seed)
  X = rmvnorm(100, rep(0, 30), block.matrix)
  
  #add colnames
  features <- c(sprintf("X%d", seq(1,30)))
  colnames(X) = features
  
  #define signal-to-noise ratio/compute variance of error
  snr = 1.8
  var.error = (t(beta) %*% block.matrix %*% beta)/snr
  
  #simulate y = X^T * beta
  set.seed(-1)
  y = as.numeric(X %*% beta + rnorm(100, 0, sqrt(var.error)))
  
  #return object
  sim.data = cbind(y,X)
  colnames(sim.data)[1] = "y"
  
  return(sim.data)
}

# create a list to contains different data based on different r values.
sim.data.list = lapply(r, function(x) sim_data(x, block.num, block.group, nrow, ncol, beta, 123))


#====#====#====#====#====#==== function for computing Sensitivity/Specificity
metrics = function(pred_vec, true_vec){
  TP = sum(pred_vec  == 1 & true_vec == 1)
  FP = sum(pred_vec  == 1 & true_vec == 0)
  FN = sum(pred_vec  == 0 & true_vec == 1)
  TN = sum(pred_vec  == 0 & true_vec == 0)
  
  sensitivity = TP/(TP+FN)
  specificity = TN/(FP+TN)
  
  return(c(sensitivity, specificity))
}

#====#====#====#====#====#==== function for best tune lambda for lasso
lambda.tune = function(x, y){
  cv.lasso <- cv.glmnet(x= x, y=y, alpha=1, type.measure = "deviance")
  return(cv.lasso$lambda.1se)
}
#====#====#====#====#====#============== Compare mean square loss for stepwise and lasso
#====#====#====#====#====#==== function for computing mse using bootstrapping
boostrap.mse = function(times, data, seed){
  #define mean square loss
  stepwise.loss = c()
  lasso.loss = c()
  
  set.seed(seed)
  for(i in 1:times){
    #define train/test data
    boost.index = sample(1:nrow(data), replace = T)
    boost.sample = as.data.frame(data[boost.index, ])
    test.sample = as.data.frame(data[-boost.index, ])
    
    #fitting stepwise and lasso
    step.model = step(lm(y ~ ., data = boost.sample), trace = F)
    
    best.tune = lambda.tune(as.matrix(boost.sample[,-1]), boost.sample$y)
    lasso.model = glmnet(x= as.matrix(boost.sample[,-1]), y = boost.sample$y, 
                         alpha = 1, lambda = best.tune)
    
    #compute the mean square loss
    stepwise.loss[i] = mean((test.sample$y - predict(step.model, test.sample))^2)
    lasso.loss[i] = mean((test.sample$y - predict(lasso.model, as.matrix(test.sample[,-1])))^2)
  }
  return(c(mean(stepwise.loss), mean(lasso.loss)))
}

# define a data frame for mean square loss
loss.data = data.frame(
  "Stepwise" = sapply(sim.data.list, function(x) boostrap.mse(100, x, 123))[1,],
  "Lasso" = sapply(sim.data.list, function(x) boostrap.mse(100, x, 123))[2,]
)
names(loss.data) = r

#====#====#====#====#====#============== variable selection with stepwise and lasso based on bootstrapping
#====#====#====#====#====#==== function for variable selection using bootstrapping
boostrap.vars= function(times, data, seed){
  #define selected vars 
  stepwise.vars = matrix(NA, rep(30*times), nrow = times, ncol = 30)
  lasso.vars = matrix(NA, rep(30*times), nrow = times, ncol = 30)
  
  #loop through each times and select vars
  set.seed(seed)
  for(i in 1:times){
    #define train/test data
    boost.index = sample(1:nrow(data), replace = T)
    boost.sample = as.data.frame(data[boost.index, ])
    
    #fitting stepwise and lasso
    step.model = step(lm(y ~ ., data = boost.sample), trace = F)
    
    best.tune = lambda.tune(as.matrix(boost.sample[,-1]), boost.sample$y)
    lasso.model = glmnet(x= as.matrix(boost.sample[,-1]), y = boost.sample$y, 
                         alpha = 1, lambda = best.tune)
    
    ####variable selection, if selected, 1, otherwise 0
    #####NOTE !!!!!
    #stepwise model only shows the active betas, but I need all the betas, 
    #if not active, it should set as 0, like lasso$beta if this confuse you.
    step.coef = coef(step(lm(y~., data = boost.sample), trace = F))[-1]
    #find the outset of total coef and step.coef
    outset = setdiff(names(true.coef.matrix[1,]), names(step.coef))
    temp = rep(c(0), c(length(outset)))
    names(temp) = outset
    #combine null set and active set together, becasue my computation needs TNR
    step.coef.tot = c(step.coef, temp)
    
    #if selected, 1, otherwise 0
    stepwise.vars[i, ] =  ifelse(step.coef.tot != 0 , 1, 0)
    lasso.vars[i, ] = ifelse(lasso.model$beta[,1] != 0 , 1, 0)
  }
  return(list(stepwise.vars, lasso.vars))
}

# define the matrix as true vars for testing 
true.coef.matrix = matrix(0, rep(30), nrow = 1, ncol = 30)
colnames(true.coef.matrix) = c(sprintf("X%d", seq(1,30)))
true.coef.matrix[, which(beta[,1] != 0)] = 1

# use lapply loop through each correlation case.
stepwise.vars.list = lapply(sim.data.list, function(x) boostrap.vars(50, x, 123)[[1]])
lasso.vars.list = lapply(sim.data.list, function(x) boostrap.vars(50, x, 123)[[2]])

# compute sensitivity/specificity
stepwise.metrics = list()
lasso.metrics = list()

for(i in 1:3){
  stepwise.metrics[[i]] = t(sapply(1:50, function(x) metrics(stepwise.vars.list[[i]][x,], true.coef.matrix[1,])))
  lasso.metrics[[i]] =  t(sapply(1:50, function(x) metrics(lasso.vars.list[[i]][x,], true.coef.matrix[1,])))
}

# compute the avg sensitivity/specificity
stepwise.TPR = sapply(1:3, function(x) mean(stepwise.metrics[[x]][,1]))
stepwise.TNR = sapply(1:3, function(x) mean(stepwise.metrics[[x]][,2]))

metrics.data = data.frame(
  "Lasso_TPR" = sapply(1:3, function(x) mean(lasso.metrics[[x]][,1])),
  "stepwise_TPR" = sapply(1:3, function(x) mean(stepwise.metrics[[x]][,1])),
  "Lasso_TNR" = sapply(1:3, function(x) mean(lasso.metrics[[x]][,2])),
  "stepwise_TNR" = sapply(1:3, function(x) mean(stepwise.metrics[[x]][,2]))
)
names(metrics.data) = r




