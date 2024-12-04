rm(list=ls()) #Clean-up working environment

##########################################################
# THIS R SCRIPT ALLOWS REPLICATING SIMULATION SCENARIO 2 #
##########################################################

###########################################
# Set up the parameters of the simulation #
###########################################

n.iter <- 5e3   # Number of simulated datasets for each combination of parameters
n      <- 20000 # Sample size

# Simulation scenarios
parameters <- data.frame(scenario=c('2a',  #No residual confounding, but the method cannot detect that due to misspecification of the model for the exposure (X)
                                    '2b')) #Residual confounding, but the method cannot detect that due to misspecification of the model for the outcome (Y)

# The object res will store the desired results from the simulations
res <- data.frame(parameters,
                  RR_X_linear=NA,    # Rejection rate for the null hypothesis test that the coefficient for covariate 1 (C1) from a linear regression model for the exposure (X) given C1 and covariate 2 (C2)
                  RR_X_quadratic=NA, # Rejection rate for the null hypothesis test that the coefficient for C1 from a linear regression model for X given C1 and C2, where C2 is modeled as a quadratic term
                  bias_linear=NA,    # Bias of the coefficient for X from a linear regression model for the outcome (Y) given X, C1 and C2
                  RR_Y_linear=NA,    # Rejection rate for the null hypothesis test that the coefficient for C1 from a linear regression model for Y given X, C1 and C2
                  bias_quadratic=NA, # Bias of the coefficient for X from a linear regression model for Y given X, C1 and C2
                  RR_Y_quadratic=NA) # Rejection rate for the null hypothesis test that the coefficient for C1 from a linear regression model for Y given X, C1 and C2, where C2 is modeled as a quadratic term

#######################
# Run the simulations #
#######################

# Loop through the two scenarios
for(a in 1:nrow(parameters)) {
 
  scenario.a <- parameters[a,]
  
  res.a <- matrix(nrow=n.iter,ncol=6) #Matrix to store desired results from each simulated dataset under the current scenario
  
  for(b in 1:n.iter) {

    if (scenario.a=='2a') {
      
      # Simulate the dataset under scenario 2a
      U  <- sample(seq(-2,2,length.out=n))
      C1 <- 0.5*U; C1 <- C1 + rnorm(n,sd=sqrt(1-var(C1)))
      C2 <- rnorm(n)
      X  <- 0.5*C2 + 0.5*U^2; X <- X + rnorm(n,sd=sqrt(1-var(X)))
      Y  <- rnorm(n)

    } else if(scenario.a=='2b') {
      
      # Simulate the dataset under scenario 2b
      U  <- rnorm(n)
      C1 <- 0.5*U; C1 <- C1 + rnorm(n,sd=sqrt(1-var(C1)))
      C2 <- rnorm(n)
      X  <- 0.5*C2 + 0.5*(U+0.2)^2; X <- X + rnorm(n,sd=sqrt(1-var(X)))
      Y  <- 0.3269946*(U+0.1417565)^2; Y <- Y + rnorm(n,sd=sqrt(1-var(Y)))
    }
    
    #Fit four linear regression models
    fitX_linear    <- summary(lm(X~C1+C2))$coef      #Model 1, linear
    fitX_quadratic <- summary(lm(X~I(C1^2)+C2))$coef #Model 1, quadratic
    
    fitY_linear    <- summary(lm(Y~X+C1+C2))$coef      #Model 2, linear
    fitY_quadratic <- summary(lm(Y~X+I(C1^2)+C2))$coef #Model 2, quadratic
    
    res.a[b,] <- c(fitX_linear[2,4],    #Extract the P-value for the coefficient for C1 from model 1, linear
                   fitX_quadratic[2,4], #Extract the P-value for the coefficient for C1 from model 1, quadratic
                   fitY_linear[2,1],    #Extract the regression coefficient for X from model 2, linear
                   fitY_linear[3,4],    #Extract the P-value for the coefficient for C1 from model 2, linear
                   fitY_quadratic[2,1], #Extract the regression coefficient for X from model 2, quadratic
                   fitY_quadratic[3,4]) #Extract the P-value for the coefficient for C1 from model 2, quadratic
    
  }

  # Summarize the results for the current scenario
  res$RR_X_linear[a]    <- mean(res.a[,1]<0.05) # Proportion of simulated datasets where the P<5% for the coefficient of C1 from model 1, linear 
  res$RR_X_quadratic[a] <- mean(res.a[,2]<0.05) # Proportion of simulated datasets where the P<5% for the coefficient of C1 from model 1, quadratic
  res$bias_linear[a]    <- mean(res.a[,3])      # Average (over simulated datasets0 of the coefficient for X from model 2, linear 
  res$RR_Y_linear[a]    <- mean(res.a[,4]<0.05) # Proportion of simulated datasets where the P<5% for the coefficient of C1 from model 2, quadratic
  res$bias_quadratic[a] <- mean(res.a[,5])      # Average (over simulated datasets0 of the coefficient for X from model 2, quadratic
  res$RR_Y_quadratic[a] <- mean(res.a[,6]<0.05) # Proportion of simulated datasets where the P<5% for the coefficient of C1 from model 2, quadratic
}

#Save the simulation results in a .txt file
write.table(res, 'Simulation_1_results.txt', row.names=F, sep='\t')
