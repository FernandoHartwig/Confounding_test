##########################################################
# THIS R SCRIPT ALLOWS REPLICATING SIMULATION SCENARIO 3 #
##########################################################

rm(list=ls()) #Clean-up working environment

#Set up parallel computing
library(doParallel)
registerDoParallel(cores=n.cores) # Replace n.cores with the number of cores too be used in the simulations

###########################################
# Set up the parameters of the simulation #
###########################################

n.iter    <- 5000 #number of simulated datasets for each combination of parameters
ncov      <- 20   #number of covariates (constant over all scenarios)
ncov_conf <- 12   #number of covariates that have a direct effect on both the exposure X and the outcome Y (constant over all scenarios)

#n:            sample size
#stat_cutoff:  statistical cutoff to classify a covariate as having satisfied conditions (i) and (ii) of the Theorem. 1: more relaxed; 2: more stringent.
#c_X:          number of covariates that have a direct effect on X, but not on Y
#c_Y:          number of covariates that have a direct effect on Y, but not on X
#rho:          Pearson correlation between each pair of covariates
#f_X:          True data-generating model for X
#f_Y:          True data-generating model for Y
#model_spec:   Assumed models for X and Y
#theta_X:      %Variance in X explained by confounders in a correctly-specified model (also regulates the effect of X-only causes on X)
#theta_Y:      %Variance in Y explained by confounders in a correctly-specified model (also regulates the effect of Y-only causes on Y)

#The parameters object will be a dataframe containing all desired combinations of parameters

#3a: Varying model specification and C_X
parameters <- data.frame(scenario='3a', n=5000, stat_cutoff=2, rho=0, c_X=rep(c(0,4),3), c_Y=4, theta_X=0.3, theta_Y=0.3, f_X=1, f_Y=1, model_spec=rep(c(1:3),each=2))

#3b: Varying stat_cutoff
parameters <- rbind(parameters,
                    data.frame(scenario='3b', n=5000, stat_cutoff=1, rho=0, c_X=4, c_Y=4, theta_X=0.3, theta_Y=0.3, f_X=1, f_Y=1, model_spec=c(1,3)))

#3c: Varying the true data-generating models
parameters <- rbind(parameters,
                    data.frame(scenario='3c', n=5000, stat_cutoff=2, rho=0, c_X=4, c_Y=4, theta_X=0.3, theta_Y=0.3, f_X=rep(c(1,2,3),each=2), f_Y=rep(c(2,3),3), model_spec=3))

#3d: Varying theta_X and theta_Y
parameters <- rbind(parameters,
                    data.frame(scenario='3d', n=5000, stat_cutoff=2, rho=0, c_X=4, c_Y=4, theta_X=rep(c(0.2,0.6),each=2), theta_Y=rep(c(0.2,0.6),2), f_X=1, f_Y=1, model_spec=3))

#3e: Varying rho
parameters <- rbind(parameters,
                    data.frame(scenario='3e', n=5000, stat_cutoff=2, rho=c(0.3,0.6), c_X=4, c_Y=4, theta_X=0.3, theta_Y=0.3, f_X=1, f_Y=1, model_spec=3))

#3f: Setting C_Y=0 (unlike other scenarios where C_Y=4)
parameters <- rbind(parameters,
                    data.frame(scenario='3f', n=5000, stat_cutoff=2, rho=0, c_X=4, c_Y=0, theta_X=0.3, theta_Y=0.3, f_X=1, f_Y=1, model_spec=3))

####################################################################################################
# If the simulations stop for some reason, one can pick up from the last combination of parameters #
####################################################################################################

filename <- 'Simulation_3_results.txt'

if(!file.exists(filename)) {

  #The object res will store the desired results from the simulations
  res <- data.frame(parameters,
                    XY_bias=NA,        #Estimated bias in the coefficient for X in a linear regression model for Y given X and covariates, averaged over all simulated datasets
                    ncovs_i=NA,        #Number of covariates satisfying condition (i) of the Theorem, averaged over all simulated datasets
                    propcovs_i_ii0=NA, #Proportion of simulated datasets where no covariate was classified as satisfying both condition (i) and condition (ii) of the Theorem
                    propcovs_i_ii1=NA, #Proportion of simulated datasets where at least one covariate was classified as satisfying both condition (i) and condition (ii) of the Theorem
                    propcovs_i_ii2=NA, #Proportion of simulated datasets where at least two covariates were classified as satisfying both condition (i) and condition (ii) of the Theorem
                    propcovs_i_ii3=NA) #Proportion of simulated datasets where at least two covariates were classified as satisfying both condition (i) and condition (ii) of the Theorem
  
  par.start <- 1
  
} else {
  
  res       <- read.table(filename, header=T, sep='\t')
  par.start <- min(which(is.na(res$XY_bias)))
}

###########################################
# Loop through combinations of parameters #
###########################################

for(a in par.start:nrow(parameters)) {
  
  #Extract current simulation parameters
  scenario.a <- parameters[a,]

  n.a            <- scenario.a$n
  stat_cutoff.a  <- scenario.a$stat_cutoff
  c_X.a          <- scenario.a$c_X
  c_Y.a          <- scenario.a$c_Y
  rho.a          <- scenario.a$rho
  f_X.a          <- scenario.a$f_X
  f_Y.a          <- scenario.a$f_Y
  model_spec.a   <- scenario.a$model_spec
  theta_X.a      <- scenario.a$theta_X
  theta_Y.a      <- scenario.a$theta_Y
  
  ######################################################################################
  # Set up objects to be used in all iterations of the a-th combinations of parameters #
  ######################################################################################
  
  if(stat_cutoff.a==1) {
    pX.a <- 0.05 #Statistical cutoff (before Bonferroni correction) to classify covariates as satisfying condition (i) of the Theorem.
    pY.a <- 0.05 #Statistical cutoff (no Bonferroni correction applied) to classify covariates as satisfying condition (ii) of the Theorem.
    
  } else {
    pX.a <- 0.001
    pY.a <- 0.3
  }
  
  #Covariates will be sampled from a multivariate normal distribution with mean vector cov_mu.a and covariance matrix cov_Sigma.a
  #Afterwards, half of the covariates will be dichotomized
  cov_mu.a          <- rep(0,ncov)
  cov_Sigma.a       <- matrix(rho.a, nrow=ncov, ncol=ncov)
  diag(cov_Sigma.a) <- 1
  
  #Set up indicators that will be shuffled in each iteration
  I_star.a    <- c(rep(1,ncov_conf),rep(0,ncov-ncov_conf))         #This will be shuffled to determine which covariates affect both exposure and outcome 
  I_X_prime.a <- c(rep(1,c_X.a),rep(0,ncov-ncov_conf-c_X.a))       #This will be shuffled to determine which covariates affect the exposure but not the outcome
  I_Y_prime.a <- c(rep(1,c_Y.a),rep(0,ncov-ncov_conf-c_X.a-c_Y.a)) #This will be shuffled to determine which covariates affect the outcome but not the exposure
  
  C_continuous_names <- paste('C_continuous', 1:(ncov/2), sep='') #Vector containing the names of continuous covariates
  C_binary_names     <- paste('C_binary', 1:(ncov/2), sep='')     #Vector containing the names of binary covariates
  C_names            <- c(C_continuous_names, C_binary_names)

  #The res.a object will store results for each dataset simulated using the a-th combination of parameters
  res.a <- foreach(b=1:n.iter, .combine='rbind', .inorder=F, .packages=c('MASS', 'sandwich', 'lmtest')) %dopar% {
    
    #######################
    # Simulate covariates #
    #######################
    
    #Generate jointly normally distributed covariates
    C <- mvrnorm(n.a, cov_mu.a, cov_Sigma.a)
    
    #Dichotomize half of the covariates
    C_continuous <- C[,1:(ncov/2)]
    C_binary     <- C[,(ncov/2+1):ncov]
    for(c in 1:ncol(C_binary)) {
      C_binary[,c] <- as.numeric(C_binary[,c]>=0)
    }
    C           <- data.frame(C_continuous, C_binary)
    colnames(C) <- C_names
    
    #Determine which covariates affect exposure, outcome or both
    I_star    <- sample(I_star.a)
    I_X_prime <- rep(0,ncov); I_X_prime[I_star==0] <- sample(I_X_prime.a)
    I_Y_prime <- rep(0,ncov); I_Y_prime[I_star==0 & I_X_prime==0] <- sample(I_Y_prime.a)
    
    #####################
    # Simulate exposure #
    #####################
    
    alpha_X <- runif(ncov/2,0.01,0.03) #Coefficients for the "main effects" of continuous covariates on X
    beta_X  <- 8*alpha_X               #Coefficients for the "main effects" of binary covariates on X
    gamma_X <- -2*beta_X               #Coefficients for product terms on X

    #Set to zero the coefficients associated with covariates that are no X-Y confounders
    alpha_X_star <- alpha_X*I_star[1:(ncov/2)]
    beta_X_star  <- beta_X*I_star[(ncov/2+1):ncov]
    
    #Set to zero the coefficients associated with covariates that affect X, but not Y
    gamma_X_star <- gamma_X*I_star[1:(ncov/2)]
    alpha_X_prime <- alpha_X*I_X_prime[1:(ncov/2)]
    beta_X_prime  <- beta_X*I_X_prime[(ncov/2+1):ncov]
    gamma_X_prime <- gamma_X*I_X_prime[1:(ncov/2)]

    X_star_tilde  <- (C_continuous^3)%*%alpha_X_star  + C_binary%*%beta_X_star  + (C_binary*C_continuous)%*%gamma_X_star  #Component of X that is influenced by X-Y confounders
    X_prime_tilde <- (C_continuous^3)%*%alpha_X_prime + C_binary%*%beta_X_prime + (C_binary*C_continuous)%*%gamma_X_prime #Component of X that is influenced by covariates that affect X, but not Y
    
    #Apply a transformation depending the the true model for X (f_X)
    if(f_X.a==1) {
      X_star  <- X_star_tilde
      X_prime <- X_prime_tilde
      
    } else if (f_X.a==2) {
      X_star  <- exp(X_star_tilde)
      X_prime <- exp(X_prime_tilde)
      
    } else if (f_X.a==3) {
      X_star  <- log(X_star_tilde-min(X_star_tilde)+2)
      X_prime <- log(X_prime_tilde-min(X_prime_tilde)+2)
    }

    X <- sqrt(theta_X.a)*X_star/sd(X_star); if(c_X.a>0) { X <- X + sqrt(c_X.a/ncov_conf*theta_X.a)*X_prime/sd(X_prime) }
    X <- X + rnorm(n.a,sd=sqrt(1-var(X))) #Exposure has unit variance

    ####################
    # Simulate outcome #
    ####################

    alpha_Y <- runif(ncov/2,0.01,0.03) #Coefficients for the "main effects" of continuous covariates on Y
    beta_Y  <- 8*alpha_Y               #Coefficients for the "main effects" of binary covariates on Y
    gamma_Y <- -2*beta_Y               #Coefficients for product terms on Y

    #Set to zero the coefficients associated with covariates that are no X-Y confounders
    alpha_Y_star <- alpha_Y*I_star[1:(ncov/2)]
    beta_Y_star  <- beta_Y*I_star[(ncov/2+1):ncov]
    gamma_Y_star <- gamma_Y*I_star[1:(ncov/2)]

    #Set to zero the coefficients associated with covariates that affect X, but not Y
    alpha_Y_prime <- alpha_X*I_Y_prime[1:(ncov/2)]
    beta_Y_prime  <- beta_X*I_Y_prime[(ncov/2+1):ncov]
    gamma_Y_prime <- gamma_X*I_Y_prime[1:(ncov/2)]
    
    Y_star_tilde  <- (C_continuous^3)%*%alpha_Y_star  + C_binary%*%beta_Y_star  + (C_binary*C_continuous)%*%gamma_Y_star  #Component of Y that is influenced by X-Y confounders
    Y_prime_tilde <- (C_continuous^3)%*%alpha_Y_prime + C_binary%*%beta_Y_prime + (C_binary*C_continuous)%*%gamma_Y_prime #Component of Y that is influenced by covariates that affect Y, but not X
    
    #Apply a transformation depending the the true model for Y (f_Y)
    if(f_Y.a==1) {
      Y_star  <- X_star_tilde
      Y_prime <- Y_prime_tilde
      
    } else if(f_Y.a==2) {
      Y_star  <- exp(X_star_tilde)
      Y_prime <- exp(Y_prime_tilde)
      
    } else if(f_Y.a==3) {
      Y_star  <- log(X_star_tilde)
      Y_prime <- log(Y_prime_tilde)
    }

    Y <- sqrt(theta_Y.a)*Y_star/sd(Y_star); if(c_Y.a>0) { Y <- Y + sqrt(c_Y.a/ncov_conf*theta_Y.a)*Y_prime/sd(Y_prime) }
    Y <- Y + rnorm(n.a,sd=sqrt(1-var(Y))) #Outcome has unit variance

    #####################################################################################
    # Classify each covariate as satisfying each conditions (i) and (ii) of the Theorem #
    #####################################################################################

    #Specify which covariates will be included and how they will be modeled    
    if(model_spec.a==1) {
      
      #All covariates, with cubic terms for continuous covariates, terms for binary covariates and 10 product terms, each involving one continuous and one binary covariate
      terms_model <- c(paste('I(', C_continuous_names, '^3)', sep=''), C_binary_names, paste(C_continuous_names, ':', C_binary_names, sep=''))
      C_names.a   <- C_names
      
    } else if(model_spec.a==2) {
      
      #All covariates, each modelled as a linear term, with no product terms or transformations of continuous covariates
      terms_model <- c(C_continuous_names, C_binary_names)
      C_names.a   <- C_names
      
    } else if(model_spec.a==3) {

      #Same as model_spec.a==2, but excluding 6 confounders from the covariate set      
      C_names_to_remove   <- sample(C_names[I_C==1],6)
      ind_keep_continuous <- !C_continuous_names%in%C_names_to_remove
      ind_keep_binary     <- !C_binary_names%in%C_names_to_remove
      terms_model         <- c(C_continuous_names[ind_keep_continuous], C_binary_names[ind_keep_binary])
      C_names.a           <- C_names[!C_names%in%C_names_to_remove]
    }
    
    ncov.a <- length(C_names.a)

    #Fit linear regression models for exposure and outcome    
    Xformula <- paste('X~',   paste(terms_model, collapse='+'), sep='')
    Yformula <- paste('Y~X+', paste(terms_model, collapse='+'), sep='')
    
    data_lm <- data.frame(Y,X,C)
    
    fitX <- lm(Xformula, data=data_lm)
    fitY <- lm(Yformula, data=data_lm)

    #For each covariate, obtain P-values to classify the covariate as satisfying condition (i) and condition (ii) of the Theorem
    if(model_spec.a==1) {

      ind_covs_i  <- ind_covs_ii <- numeric(ncov.a)
      
      vcovHC_X <- vcovHC(fitX)
      vcovHC_Y <- vcovHC(fitY)
      
      for(j in 1:ncov.a) {
        
        cov_type   <- c('continuous','binary')[c(grepl('continuous',C_names.a[j]),grepl('binary',C_names.a[j]))]
        cov_number <- strsplit(C_names.a[j], cov_type)[[1]][2]
        
        coefs_to_remove <- c(paste('I(',C_names.a[j],'^3)',sep=''), paste('C_binary', cov_number, ':', C_names.a[j], sep=''),
                             C_names.a[j], paste(C_names.a[j], ':', 'C_continuous', cov_number, sep=''))
        
        coefs_to_remove <- paste(coefs_to_remove, sep='', collapse='-')
        
        fitX.reduced <- update(fitX, paste('~.-', coefs_to_remove, sep=''))
        fitY.reduced <- update(fitY, paste('~.-', coefs_to_remove, sep=''))
        
        ind_covs_i[j]  <- waldtest(fitX.reduced, fitX, vcov=vcovHC_X)[2,4]<(pX.a/ncov.a)
        ind_covs_ii[j] <- waldtest(fitY.reduced, fitY, vcov=vcovHC_Y)[2,4]>=pY.a
      }
      
    } else {
      ind_covs_i  <- coeftest(fitX, vcov=vcovHC(fitX))[-1,4]<(pX.a/ncov.a)
      ind_covs_ii <- coeftest(fitY, vcov=vcovHC(fitY))[-c(1,2),4]>=pY.a
    }
    
    ############################################
    # Extract results from the current dataset #
    ############################################
    
    XY_bias       <- fitY$coefficients[2]
    ncovs_i       <- sum(ind_covs_i)
    ncovs_i_ii    <- sum((ind_covs_i+ind_covs_ii)==2)
    
    c(XY_bias, ncovs_i, ncovs_i_ii)
  }
  
  ############################################################
  # Summarize results for the a-th combination of parameters #
  ############################################################
  
  res$XY_bias[a]                <- mean(res.a[,1])
  res$ncovs_i[a]                <- mean(res.a[,2])
  res$propcovs_i_ii0[a]         <- mean(res.a[,3]==0)
  res$propcovs_i_ii1[a]         <- mean(res.a[,3]>=1)
  res$propcovs_i_ii2[a]         <- mean(res.a[,3]>=2)
  res$propcovs_i_ii3[a]         <- mean(res.a[,3]>=3)

  print(paste('Progress: ', round(100*a/nrow(parameters),2), '%.', sep=''))
  
  write.table(res, filename, row.names=F, sep='\t')
}