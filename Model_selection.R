####################################################################################################################################################
# THIS R SCRIPT ALLOWS REPLICATING THE PROCESS USED TO SELECT THE SPECIFICATION FOR THE MODELS FOR EXPOSURE X AND OUTCOME Y IN THE APPLIED EXAMPLE #
####################################################################################################################################################

rm(list=ls()) #Clean-up working environment

#Load required packages
library(gtools)

##############################################
# INSTRUCTIONS ABOUT THE DATA TO BE ANALYZED #
##############################################

#Import the data into a dataframe named real.data.
#Create a vector covariate_names with the names of the covariates.
#Ensure the exposure is labelled as "X" (without the quotation marks) in the real.data data.frame, and is a binary variable coded as 0 or 1.
#Ensure the outcome is labelled as "Y" (without the quotation marks) in the real.data data.frame, and is a continuous variable.
#Ensure there is no missing data for X, Y or any covariate.

ncov <- length(covariate_names) #Number of covariates

######################################################
# SET UP FORMULAS FOR DIFFERENT MODEL SPECIFICATIONS #
######################################################

#Store all possible second-order terms involving covariates
second_order <- NULL
for(a in 1:(length(covariate_names)-1)) {
  
  if(a!=length(covariate_names)) {
    second_order <- c(second_order, paste(covariate_names[a],covariate_names[(a+1):length(covariate_names)],sep=':'))
  }
  
  if(class(real.data[,covariate_names[a]])%in%c('numeric','integer')) {
    second_order <- c(second_order, paste('I(',covariate_names[a], '^2)', sep=''))
  }
}

#Formulas describing the specification of logistic regression models for X
formulaX_list      <- list(paste('X',paste(covariate_names,collapse='+'),sep='~')) #Simple specification, with no second-order terms (1.1 in Supplementary Section 7)
formulaX_list[[2]] <- paste(c(formulaX_list[[1]], second_order),collapse='+')      #Specification with first and second-order terms (1.2 in Supplementary Section 7)
formulaX_list[[3]] <- 'X~cov_catX_aic'                                             #Specification including a single categorical variable (catX_aic) generated below by conditional inference tree analyses using the Akaike information criterion (AIC) as a stopping criterion (1.3 in Supplementary Section 7)
formulaX_list[[4]] <- 'X~cov_catX_bic'                                             #Same as above, but using the Bayesian information criterion (BIC) as a stopping criterion (1.4 in Supplementary Section 7)
formulaX_list[[5]] <- paste(c(formulaX_list[[1]], 'cov_catX_aic'),collapse='+')    #Same as 1.1, but also including cov_catX_aic
formulaX_list[[6]] <- paste(c(formulaX_list[[1]], 'cov_catX_bic'),collapse='+')    #Same as 1.1, but also including cov_catX_bic

#Formulas describing the specification of linear regression models for Y
formulaY_list      <- list(paste(c('X',covariate_names),collapse='+'))                              #Simple specification, with no second-order terms (2.1 in Supplementary Section 7)
formulaY_list[[2]] <- paste(c(formulaY_list[[1]], second_order),collapse='+')                       #Specification with first and second-order terms for covariates, and only a first order term for the exposure (2.2 in Supplementary Section 7)
formulaY_list[[3]] <- paste(c(formulaY_list[[1]], paste('X',covariate_names,sep=':')),collapse='+') #Same as 2.2, but including product terms between the exposure and each covariate (2.3 in Supplementary Section 7)
formulaY_list[[4]] <- 'X+cov_catY_aic'                                                              #Specification including the exposure and a single categorical variable (catY_aic) generated below by conditional inference tree analyses using the  AIC as a stopping criterion (2.4 in Supplementary Section 7)
formulaY_list[[5]] <- 'X*cov_catY_aic'                                                              #Same as 2.4, but including a product term (2.5 in Supplementary Section 7)
formulaY_list[[6]] <- 'X+cov_catY_bic'                                                              #Same as 2.4, but usin the BIC instead of the AIC (2.6 in Supplementary Section 7)
formulaY_list[[7]] <- 'X*cov_catY_bic'                                                              #Same as 2.5, but usin the BIC instead of the AIC (2.7 in Supplementary Section 7)

#######################################################
# PERFORM MODDEL SELECTION BY 5-FOLD CROSS-VALIDATION #
#######################################################

res_cv_x <- res_cv_y <- NULL #res_cv_x will store the resuls for the models for X; res_cv_y will store the resuls for the models for Y

#For greater stability, the cross-validation process is repeated 50 times and the results are averaged over the 50 iterations
count <- 0
while(count<50) {

  count <- count+1
  
  print(count)
  
  #Define the partition of the real.data for the current cross-validation iteration
  ind_cv <- sample(as.numeric(quantcut(1:nrow(real.data),5)))

  #####################################
  # ASSESS THE DIFFERENT MODELS FOR X #
  #####################################
  
  res_cv <- NULL #This will store the results for the models for X for the current cross-validation iteration
  
  #Loop over all subsets of the real.data
  for(i in unique(ind_cv)) {
    
    #Partition the real.data into training and testing subsets
    training.i <- ind_cv!=i
    test.i     <- ind_cv==i
    
    real.data.training.i <- real.data[training.i,]
    real.data.test.i     <- real.data[test.i,]

    #Conditional inference tree analyses  
    p          <- 0.9                                                                                                                     #Starting p-value threshold (the larger the p-value, the greater the tendency to partition into more nodes - better fit and higher risk of overfitting)
    formula_X  <- formula(paste('X', paste(covariate_names, collapse='+'), sep='~'))                                                      #Formula indicating the covariates to be used
    minbucket  <- ceiling(nrow(real.data.training.i)*0.01)                                                                                #Avoid having terminal nodes with a very small number of individuals
    cur_tree   <- ctree(formula_X, data=real.data.training.i, control=ctree_control(alpha=p, minbucket=minbucket, testtype='Univariate')) #Fit the most flexible tree
    cur_fitX   <- glm(X~factor(predict(cur_tree, type='node')), data=real.data.training.i, family=binomial)                               #Fit a model for X given node membership
    AIC_values <- AIC(cur_fitX)                                                                                                           #Measure the quality of fit using the Akaike Information Criterion (AIC)

    #Loop over decreasing p-values (i.e., increasingly stringent partitioning criterion) until the AIC starts to increase
    continue <- T
    while(continue) {
      
      if(p>0.01) {
        p <- p-0.01
        
      } else {
        p <- p/2
      }
      
      next_tree  <- ctree(formula_X, data=real.data.training.i, control=ctree_control(alpha=p, minbucket=minbucket, testtype='Univariate')) #Fit a tree using the current p-value threshold under consideration
      next_fitX  <- glm(X~factor(predict(next_tree, type='node')), data=real.data.training.i, family=binomial)                              #Fit a model for X given node membership
      AIC_values <- c(AIC_values, AIC(next_fitX))                                                                                           #Measure the quality of fit using the Akaike Information Criterion (AIC)

      #Check if the AIC went up in comparison to the previous tree. If so, stop the process. If not, continue the process.
      if(AIC_values[length(AIC_values)]>AIC_values[length(AIC_values)-1]) {
        continue <- F
        
      } else {
        cur_tree <- next_tree
        cur_fitX <- next_fitX
        
      }
    }

    cur_tree_aic <- cur_tree #Store the best tree according to the AIC criterion (i.e., the three with the smallest AIC)

    #Now, continue the process, but this time using the Bayesian Information Criterion (BIC) as the stopping criterion.  
    BIC_values <- BIC(cur_fitX)
    
    continue <- T
    while(continue) {
      
      if(p>0.01) {
        p <- p-0.01
        
      } else {
        p <- p/2
      }
      
      next_tree  <- ctree(formula_X, data=real.data.training.i, control=ctree_control(alpha=p, minbucket=minbucket, testtype='Univariate'))
      next_fitX  <- glm(X~factor(predict(next_tree, type='node')), data=real.data.training.i, family=binomial)
      BIC_values <- c(BIC_values, BIC(next_fitX))
      
      if(BIC_values[length(BIC_values)]>BIC_values[length(BIC_values)-1]) {
        continue <- F
        
      } else {
        cur_tree <- next_tree
        cur_fitX <- next_fitX
        
      }
    }
    
    cur_tree_bic <- cur_tree #Store the best tree according to the BIC criterion (i.e., the three with the smallest BIC)

    #Save the AIC and BIC-optimizing trees generated using the training data
    cov_catX_aic_training <- factor(predict(cur_tree_aic, type='node'))
    cov_catX_bic_training <- factor(predict(cur_tree_bic, type='node'))
    cov_catX_training     <- data.frame(cov_catX_aic=cov_catX_aic_training, cov_catX_bic=cov_catX_bic_training)
    
    #Apply the partition rules obtained from the AIC- and BIC-optimizing trees to the test data to generate nodes using the same partition rules
    cov_catX_aic_test <- factor(predict(cur_tree_aic, real.data.test.i, type='node'))
    cov_catX_bic_test <- factor(predict(cur_tree_bic, real.data.test.i, type='node'))
    cov_catX_test     <- data.frame(cov_catX_aic=cov_catX_aic_test, cov_catX_bic=cov_catX_bic_test)

    res.i <- NULL #Define object to store the log-likelihood quantities
    
    #Loop through all model specifications for X and, for each one, calculate both the internal (i.e., in the training dataset) and external (i.e., in the test dataset) log-likelihood of the model
    for(j in 1:length(formulaX_list)) {
      fitX_j <- glm(formulaX_list[[j]], family=binomial, data=data.frame(cov_catX_training, real.data.training.i)) #Internal log-likehood
      X_j    <- predict(fitX_j, data.frame(cov_catX_test, real.data.test.i))                                       #Calculate predicted values of logit(X) by applying the training model to the test data
      res.i  <- c(res.i, logLik(glm(real.data.test.i$X~X_j, family=binomial)))                                     #External log-likehood
    }
    
    res_cv <- rbind(res_cv,res.i) #Store the log-likelihood estimates obtained in the current iteration of the cross-validation process
  }
  
  res_cv_x <- rbind(res_cv_x,
                    apply(res_cv,2,mean)) #Average over the five external log-likehood estimates to obtain a final estimate of the external validity of each model specification for X in the current iteration of the cross-validation process
  
  #####################################
  # ASSESS THE DIFFERENT MODELS FOR Y #
  #####################################
  
  res_cv <- NULL #This will store the results for the models for Y for the current cross-validation iteration
  
  #Loop over all subsets of the real.data
  for(i in unique(ind_cv)) {
    
    #Partition the real.data into training and testing subsets
    training.i <- ind_cv!=i
    test.i     <- ind_cv==i
    
    real.data.training.i <- real.data[training.i,]
    real.data.test.i     <- real.data[test.i,]

    #Conditional inference tree analyses  
    p          <- 0.9                                                                                                                     #Starting p-value threshold (the larger the p-value, the greater the tendency to partition into more nodes - better fit and higher risk of overfitting)
    formula_X  <- formula(paste('X', paste(covariate_names, collapse='+'), sep='~'))                                                      #Formula indicating the covariates to be used
    minbucket  <- ceiling(nrow(real.data.training.i)*0.01)                                                                                #Avoid having terminal nodes with a very small number of individuals
    cur_tree   <- ctree(formula_X, data=real.data.training.i, control=ctree_control(alpha=p, minbucket=minbucket, testtype='Univariate')) #Fit the most flexible tree
    cur_fitX   <- glm(X~factor(predict(cur_tree, type='node')), data=real.data.training.i, family=binomial)                               #Fit a model for X given node membership
    AIC_values <- AIC(cur_fitX)                                                                                                           #Measure the quality of fit using the Akaike Information Criterion (AIC)
    
    #Conditional inference tree analyses
    p <- 0.9                                                                                                                              #Starting p-value threshold (the larger the p-value, the greater the tendency to partition into more nodes - better fit and higher risk of overfitting)
    formula_Y  <- formula(paste('Y', paste(covariate_names, collapse='+'), sep='~'))                                                      #Formula indicating the covariates to be used
    minbucket  <- ceiling(nrow(real.data.training.i)*0.01)                                                                                #Avoid having terminal nodes with a very small number of individuals
    cur_tree   <- ctree(formula_Y, data=real.data.training.i, control=ctree_control(alpha=p, minbucket=minbucket, testtype='Univariate')) #Fit the most flexible tree
    cur_fitY   <- lm(Y~factor(predict(cur_tree, type='node')), data=real.data.training.i)                                                 #Fit a model for Y given node membership
    AIC_values <- AIC(cur_fitY)                                                                                                           #Measure the quality of fit using the AIC

    #Loop over decreasing p-values (i.e., increasingly stringent partitioning criterion) until the AIC starts to increase
    continue <- T
    while(continue) {
      
      if(p>0.01) {
        p <- p-0.01
        
      } else {
        p <- p/2
      }
      
      next_tree  <- ctree(formula_Y, data=real.data.training.i, control=ctree_control(alpha=p, minbucket=minbucket, testtype='Univariate')) #Fit a tree using the current p-value threshold under consideration
      next_fitY  <- lm(Y~factor(predict(next_tree, type='node')), data=real.data.training.i)                                                #Fit a model for Y given node membership
      AIC_values <- c(AIC_values, AIC(next_fitY))                                                                                           #Measure the quality of fit using the Akaike Information Criterion (AIC)
      
      #Check if the AIC went up in comparison to the previous tree. If so, stop the process. If not, continue the process.
      if(AIC_values[length(AIC_values)]>AIC_values[length(AIC_values)-1]) {
        continue <- F
        
      } else {
        cur_tree <- next_tree
        cur_fitY <- next_fitY
      }
    }

    cur_tree_aic <- cur_tree #Store the best tree according to the AIC criterion (i.e., the three with the smallest AIC)
    
    #Now, continue the process, but this time using the Bayesian Information Criterion (BIC) as the stopping criterion.  
    BIC_values <- BIC(cur_fitY)
    
    continue <- T
    while(continue) {
      
      if(p>0.01) {
        p <- p-0.01
        
      } else {
        p <- p/2
      }
      
      next_tree  <- ctree(formula_Y, data=real.data.training.i, control=ctree_control(alpha=p, minbucket=minbucket, testtype='Univariate'))
      next_fitY  <- lm(Y~factor(predict(next_tree, type='node')), data=real.data.training.i)
      BIC_values <- c(BIC_values, BIC(next_fitY))
      
      if(BIC_values[length(BIC_values)]>BIC_values[length(BIC_values)-1]) {
        continue <- F
        
      } else {
        cur_tree <- next_tree
        cur_fitY <- next_fitY
      }
    }

    cur_tree_bic <- cur_tree #Store the best tree according to the BIC criterion (i.e., the three with the smallest BIC)

    #Save the AIC and BIC-optimizing trees generated using the training data
    cov_catY_aic      <- factor(predict(cur_tree_aic, type='node'))
    cov_catY_bic      <- factor(predict(cur_tree_aic, type='node'))
    cov_catX_training <- data.frame(cov_catX_aic=cov_catX_aic_training, cov_catX_bic=cov_catX_bic_training)
    
    #Apply the partition rules obtained from the AIC- and BIC-optimizing trees to the test data to generate nodes using the same partition rules
    cov_catY_aic_test <- factor(predict(cur_tree_aic, real.data.test.i, type='node'))
    cov_catY_bic_test <- factor(predict(cur_tree_aic, real.data.test.i, type='node'))
    cov_catY_test     <- data.frame(cov_catY_aic=cov_catY_aic_test, cov_catY_bic=cov_catY_bic_test)
    
    res.i <- NULL #Define object to store the log-likelihood quantities
    
    #Loop through all model specifications for X and, for each one, calculate both the internal (i.e., in the training dataset) and external (i.e., in the test dataset) log-likelihood of the model
    for(j in 1:length(formulaY_list)) {
      fitY_j <- lm(paste('Y', formulaY_list[[j]], sep='~'), data=daa.frame(cov_catY_training, real.data.training.i)) #Internal log-likehood
      Y_j    <- predict(fitY_j, data.frame(cov_catY_test,real.data.test.i))                                          #Calculate predicted values of Y by applying the training model to the test data
      res.i  <- c(res.i, logLik(lm(real.data.test.i$Y~Y_j)))                                                         #External log-likehood
    }
    
    res_cv <- rbind(res_cv,res.i) #Store the log-likelihood estimates obtained in the current iteration of the cross-validation process
    
  }
  
  res_cv_y <- rbind(res_cv_y,
                    apply(res_cv,2,mean)) #Average over the five external log-likehood estimates to obtain a final estimate of the external validity of each model specification for X in the current iteration of the cross-validation process
  
}

#Average over all iterations to obtain final estimate of th external validity of each model specification for X and Y
res_cv_x <- round(c(apply(res_cv_x,2,mean), NA),1)
res_cv_y <- round(apply(res_cv_y,2,mean),1)

model_x <- c(paste('1.', 1:length(formulaX_list), sep=''), NA) #Vector of the numbering scheme of the models for X
model_y <- paste('2.', 1:length(formulaY_list), sep='')        #Vector of the numbering scheme of the models for Y

res_cv <- data.frame(model_x,res_cv_x,model_y,res_cv_y) #Dataframe containing cross-validation results

#Save the cross-validation results in .txt files
write.table(res_cv, 'CV_results.txt', row.names=F, sep='\t')
