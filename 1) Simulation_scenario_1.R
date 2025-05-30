##########################################################
# THIS R SCRIPT ALLOWS REPLICATING SIMULATION SCENARIO 1 #
##########################################################

rm(list=ls()) #Clean-up working environment

###########################################
# Set up the parameters of the simulation #
###########################################

n.iter    <- 5e3                      # number of simulated datasets for each combination of parameters
n         <- c(500,1000,5000,20000)   # sample size
beta_c1c2 <- 0.3                      # beta_c1c2: direct effect of measured covariate 1 (C1) on measured covariate 2 (C2)
beta_c1x  <- c(0.1,0.2,0.3,0.4,0.5)   # beta_c1x:  direct effect of C1 on the exposure (X)
beta_c2x  <- 0.3                      # direct effect of C2 on X
beta_c2y  <- 0.3                      # direct effect of C2 on the outcome (Y)
beta_u    <- c(0,0.1,0.2,0.3,0.4,0.5) # direct effect of unmeasured variable (U) on both X and Y

# The parameters object will be a dataframe containing all desired combinations of parameters

# Scenario 1a: vary n, keeping all other parameters constant
parameters <- data.frame(scenario='1a', n, beta_c1c2, beta_c1x=0.3, beta_c2x, beta_c2y, beta_u=0.3)

# Scenario 2a: vary beta_c1x, keeping all other parameters constant
parameters <- rbind(parameters, data.frame(scenario='1b', n=5000, beta_c1c2, beta_c1x, beta_c2x, beta_c2y, beta_u=0.3))

# Scenario 3a: vary beta_u, keeping all other parameters constant
parameters <- rbind(parameters, data.frame(scenario='1c', n=5000, beta_c1c2=0.3, beta_c1x=0.3, beta_c2x, beta_c2y, beta_u))

# The object res will store the desired results from the simulations
res <- data.frame(parameters,
                  p5=NA,  # Proportion of simulated datasets where a T-test for the coefficient for C1 in a linear regression model for Y given X, C1 and C2 yielded a P-value < 5%
                  p1=NA,  # Same as p5, but using a 1% threshold instead
                  p01=NA) # Same as p1, but using a 0.1% threshold instead

#######################
# Run the simulations #
#######################

# Loop through combinations of parameters
for(a in 1:nrow(parameters)) {
 
  parameters.a <- parameters[a,]
  
  scenario.a  <- parameters.a$scenario
  n.a         <- parameters.a$n
  beta_c1c2.a <- parameters.a$beta_c1c2
  beta_c1x.a  <- parameters.a$beta_c1x
  beta_c2x.a  <- parameters.a$beta_c2x
  beta_c2y.a  <- parameters.a$beta_c2y
  beta_u.a    <- parameters.a$beta_u

  res.a <- numeric(n.iter) #Vector to store desired results from each simulated dataset under the current combination of parameters

  # Simulate n.iter datasets under the current combination of parameters
  for(b in 1:n.iter) {

    # Simulate the dataset
    C1 <- rnorm(n.a)
    U  <- rnorm(n.a)
    C2 <- beta_c1c2.a*C1; C2 <- C2 + rnorm(n.a,sd=sqrt(1-var(C2)))
    X  <- beta_c1x.a*C1 + beta_c2x.a*C2 + beta_u.a*U; X <- X + rnorm(n.a,sd=sqrt(1-var(X)))
    Y  <- beta_c2y.a*C2 + beta_u.a*U; Y <- Y + rnorm(n.a,sd=sqrt(1-var(Y)))
    
    fit <- lm(Y~X+C1+C2) #Linear regression model for Y given X, C1 and C2
    
    res.a[b,] <- summary(fit)$coef[3,4] #Extract the P-value for the coefficient for C1
  }
  
  # Summarize the results for the current combination of parameters
  res$p5[a]  <- mean(res.a<0.05)
  res$p1[a]  <- mean(res.a<0.01)
  res$p01[a] <- mean(res.a<0.001)
  
}

#Save the simulation results in a .txt file
write.table(res, 'Simulation_1_results.txt', row.names=F, sep='\t')
