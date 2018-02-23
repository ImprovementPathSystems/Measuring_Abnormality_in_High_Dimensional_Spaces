#was getting errors when using the regular when running the nBartlett function in the nFactors package in parallel. 
Custom_nBartlett <- function (x, N, alpha = 0.05, cor = TRUE,  correction = TRUE){
  
  n <- length(x)
  detail <- NULL
  bartlett.n <- anderson.n <- lawley.n <- 0
  bartlett <- bartlett.chi <- bartlett.df <- bartlett.p <- numeric(n)
  anderson.chi <- anderson.df <- anderson.p <- numeric(n)
  lawley.chi <- lawley.df <- lawley.p <- numeric(n)
  for (k in 0:(n - 1)) {
    i <- k + 1
    bartlett[i] <- prod(x[(k + 1):n])/(sum(x[(k + 1):n])/(n - 
                                                            k))^(n - k)
    bartlett.chi[i] <- -(N - 1 - ((2 * n + 5)/6) - ((2 * 
                                                       k)/3)) * log(bartlett[i])
    bartlett.df[i] <- 0.5 * (n - k) * (n - k - 1)
    if (correction == TRUE & bartlett.n > 0) 
      bartlett.df[i] <- 0.5 * (n - k + 2) * (n - k - 1)
    bartlett.p[i] <- pchisq(bartlett.chi[i], bartlett.df[i], 
                            lower.tail = FALSE)
    anderson.chi[i] <- -N * log(bartlett[i])
    anderson.df[i] <- 0.5 * (n - k + 2) * (n - k - 1)
    anderson.p[i] <- pchisq(anderson.chi[i], anderson.df[i], 
                            lower.tail = FALSE)
    lMean <- mean(x[(k + 1):n])
    lawley.chi[i] <- -(N - 1 - ((2 * n + 5)/6) - ((2 * k)/3) + 
                         sum((lMean^2)/((x[k] + lMean)^2))) * log(bartlett[i])
    lawley.df[i] <- 0.5 * (n - k) * (n - k - 1)
    lawley.p[i] <- pchisq(lawley.chi[i], lawley.df[i], lower.tail = FALSE)
    if (i == 1) {
      bartlett.n <- bartlett.n + as.numeric(bartlett.p[i] <= 
                                              alpha)
      anderson.n <- anderson.n + as.numeric(anderson.p[i] <= 
                                              alpha)
      lawley.n <- lawley.n + as.numeric(lawley.p[i] <= 
                                          alpha)
    }
    if (i > 1) {
      if (bartlett.p[i - 1] <= 0.05) 
        bartlett.n <- bartlett.n + as.numeric(bartlett.p[i] <= 
                                                alpha)
      if (anderson.p[i - 1] <= 0.05) 
        anderson.n <- anderson.n + as.numeric(anderson.p[i] <= 
                                                alpha)
      if (lawley.p[i - 1] <= 0.05) 
        lawley.n <- lawley.n + as.numeric(lawley.p[i] <= 
                                            alpha)
    }
  }
  
  if (bartlett.n == 0)  bartlett.n <- n
  # if (anderson.n == 0) anderson.n <- n
  if (lawley.n == 0) lawley.n <- n
  
  res <- list(Factors = c(bartlett = bartlett.n, 
                          anderson = anderson.n, lawley = lawley.n))
  return(res)
}

# Broken-stick function
#https://alstatr.blogspot.com/2014/12/principal-component-analysis-on-imaging.html
brStick <- function (x) {
  m <- 0
  out <- matrix(NA, ncol = 2, nrow = length(x))
  colnames(out) <- c("% of Variability", "B-Stick Threshold")
  for (i in 1:length(x)) {
    for (k in i:length(x)) {
      m <- m + ((1 / length(x)) * (1 / k))
    }
    out[i, ] <- c((x[i] / sum(x)) * 100, m * 100)
    m <- 0
  }
  return(max(which(out[, 1] > out[, 2])))
}


#this is a slighly modified GGI since the euclidean distance is standardized to the number of variables. I also include an option to simply calculate a euclidean distance of the original data instead of the PCA projections. 
Euclidean_Dist <- function(Patient_TS, Norm_TS,Min_VAF = .98, NoPCA = FALSE){
  if(NoPCA){
    #Compute the mean and std dev of each variable to use later to center and scale the subject data
    n_vars <- ncol(Patient_TS)
    refpop_means <- colMeans(Norm_TS)
    refpop_sd <- apply(Norm_TS,2,sd)
    Patient_TS <- data.frame(t(Patient_TS))
    dist <- (Patient_TS-refpop_means)/refpop_sd
    dist_squared <- (dist)^2
    return (sqrt(sum(dist_squared)))
    
  }else{
    #Compute the mean and std dev of each variable to use later to center and scale the subject data
    refpop_means <- colMeans(Norm_TS)
    refpop_sd <- apply(Norm_TS,2,sd)
    
    #Scaled and center. These are Z scores now.
    refpop_sc <-scale(Norm_TS, center = TRUE, scale = TRUE)
    
    #run the principal component analysis on the temporal spatial data.
    p <- prcomp(refpop_sc, center=F, scale=F, retx=T)
    
    #Compute covariance matrix
    refpop_cov<- cov(refpop_sc)
    
    #manually compute the eigenvectors and eigenvalues
    refpop_eig <- eigen(refpop_cov,symmetric=T)
    
    #Determine number of PC's to keep based on Min_VAF
    EigenSum <- sum(refpop_eig$values)
    EigenValues <- refpop_eig$values
    PC_PercentVarExplained <- t(matrix(abs(EigenValues/EigenSum)))
    numPCs<- 1
    repeat{
      if(sum(PC_PercentVarExplained[1:numPCs])>=Min_VAF){
        break
      }
      numPCs<- numPCs+ 1
    }
    
    #Calculate projections
    refpop_projs<- refpop_sc %*% p$rotation[,1:numPCs]
    
    # Get the means of the PCs for the reference groups of normal subjects
    # These should all be near zero.
    PCmeans<- t(as.matrix(colMeans(refpop_projs)))
    
    # Get the SDs of the Time Distance PCs. Should = sqrt(Eigenvalues)
    PCsds<- t(as.matrix(apply(refpop_projs,2,sd)))  
    
    Patient_TS <- data.frame(t(Patient_TS))
    Patient_TS <- data.frame(scale(Patient_TS, center = refpop_means, scale = refpop_sd))
    
    #Convert to a matrix for PC projection
    glpopm_sc<- data.matrix(Patient_TS)
    
    #Multiply the patient data by the eigenvectors calculated on the normal population. Result is the PC values for the given subject.
    glpopm_projs<- glpopm_sc%*% p$rotation[,1:numPCs]
    
    dist <- (glpopm_projs-PCmeans)/PCsds
    dist_squared <- (dist)^2
    return (sqrt(sum(dist_squared)))
  }
  
}


MAD <- function(Patient_TS, Norm_TS,Min_VAF = .98, NoPCA = FALSE){
  if(NoPCA){
    #Compute the mean and std dev of each variable to use later to center and scale the subject data
    n_vars <- ncol(Patient_TS)
    refpop_means <- colMeans(Norm_TS)
    refpop_sd <- apply(Norm_TS,2,sd)
    Patient_TS <- data.frame(t(Patient_TS))
    dist <- (Patient_TS-refpop_means)/refpop_sd
    abs_dist <- abs(dist)
    return (sum(abs_dist)/n_vars)
    
  }else{
    #Compute the mean and std dev of each variable to use later to center and scale the subject data
    refpop_means <- colMeans(Norm_TS)
    refpop_sd <- apply(Norm_TS,2,sd)
    
    #Scaled and center. These are Z scores now.
    refpop_sc <-scale(Norm_TS, center = TRUE, scale = TRUE)
    
    #run the principal component analysis on the temporal spatial data.
    p <- prcomp(refpop_sc, center=F, scale=F, retx=T)
    
    #Compute covariance matrix
    refpop_cov<- cov(refpop_sc)
    
    #manually compute the eigenvectors and eigenvalues
    refpop_eig <- eigen(refpop_cov,symmetric=T)
    
    #Determine number of PC's to keep based on Min_VAF
    EigenSum <- sum(refpop_eig$values)
    EigenValues <- refpop_eig$values
    PC_PercentVarExplained <- t(matrix(abs(EigenValues/EigenSum)))
    numPCs<- 1
    repeat{
      if(sum(PC_PercentVarExplained[1:numPCs])>=Min_VAF){
        break
      }
      numPCs<- numPCs+ 1
    }
    
    #Calculate projections
    refpop_projs<- refpop_sc %*% p$rotation[,1:numPCs]
    
    # Get the means of the PCs for the reference groups of normal subjects
    # These should all be near zero.
    PCmeans<- t(as.matrix(colMeans(refpop_projs)))
    
    # Get the SDs of the Time Distance PCs. Should = sqrt(Eigenvalues)
    PCsds<- t(as.matrix(apply(refpop_projs,2,sd)))  
    
    Patient_TS <- data.frame(t(Patient_TS))
    Patient_TS <- data.frame(scale(Patient_TS, center = refpop_means, scale = refpop_sd))
    
    #Convert to a matrix for PC projection
    glpopm_sc<- data.matrix(Patient_TS)
    
    #Multiply the patient data by the eigenvectors calculated on the normal population. Result is the PC values for the given subject.
    glpopm_projs<- glpopm_sc%*% p$rotation[,1:numPCs]
    
    dist <- (glpopm_projs-PCmeans)/PCsds
    abs_dist <- abs(dist)
    return (sum(abs_dist)/numPCs)
  }
  
}

#a function to create the principal component lines in the original space. 
Principal_Component_Lines <- function(df,SDs = 1, unit_length = FALSE){
  #SDs defines how much to scale the eigenvectors by the sqrt of their eigenvalues
  
  #scale and center
  df_sc <-scale(df, center = TRUE, scale = TRUE)
  
  #manually compute the eigenvectors and eigenvalues
  df_eig <- eigen(cov(df_sc),symmetric=T)
  evs <- sqrt(df_eig$values) #since we are taking the square root, these are standard deviations of the eigenvectors. 
  evecs <- df_eig$vectors
  
  ev1 <- evs[1]
  ev2 <- evs[2]
  
  if(evecs[ , 1][1] <0 & evecs[ , 1][2]<0){
    evecs[ , 1] <- evecs[ , 1]*-1
  }
  
  if(unit_length) {
    ev1 <- 1 ; ev2 <- 1
  }
  
  evec1_x <- c(SDs*ev1 * evecs[ , 1][1],SDs*-ev1 * evecs[ , 1][1])
  evec1_y <- c(SDs*ev1 * evecs[ , 1][2],SDs*-ev1 * evecs[ , 1][2])
  evec1_df <- data.frame(x=evec1_x,y=evec1_y)
  
  evec2_x <- c(SDs*ev2 * evecs[ , 2][1],SDs*-ev2 * evecs[ , 2][1])
  evec2_y <- c(SDs*ev2 * evecs[ , 2][2],SDs*-ev2 * evecs[ , 2][2])
  evec2_df <- data.frame(x=evec2_x,y=evec2_y)
  
  evec1_df
  
  Get_SD_Points <- function(i){
    x1 <- c(i*ev1 * evecs[ , 1][1],i*-ev1 * evecs[ , 1][1])
    y1 <- c(i*ev1 * evecs[ , 1][2],i*-ev1 * evecs[ , 1][2])
    eig1 <- rep(ev1,2)
    eig1_scale <- rep(i,2)
    x2 <- c(i*ev2 * evecs[ , 2][1],i*-ev2 * evecs[ , 2][1])
    y2 <- c(i*ev2 * evecs[ , 2][2],i*-ev2 * evecs[ , 2][2])
    eig2 <- rep(ev2,2)
    eig2_scale <- rep(i,2)
    data.frame(x=c(x1,x2), y=c(y1,y2),eig=c(eig1,eig2), eig_scale=c(eig1_scale,eig2_scale))
  }
  
  SD_Points <- lapply(1:SDs,Get_SD_Points)
  SD_Points <- do.call("rbind",SD_Points)
  
  
  return_df <- list()
  return_df$evec1_df <- evec1_df
  return_df$evec2_df <- evec2_df
  return_df$SD_Points <- SD_Points
  return_df$eig_vects <- evecs
  return_df
}



#a function to create an ellipse n standard deviations away from the mean of a 2 dimensional normal distribution. 
sdellipse <-  function(points, stdev = 1.96, density = .01, means = NULL){
  if (ncol (points) != 2) stop ('Points input must have exactly two columns.')
  if (!is.null(means) & nrow(points) > 2) stop ('Covariance matrix must be 2 by 2.')
  if (!is.null(means) & length(means) > 2) stop ('Exactly two means must be specified.')
  
  t = seq (0,2*pi+density,density)  
  x = rbind (cos(t), sin(t))
  if (is.null(means)){
    sigma = var (points)
  }
  
  A = eigen(sigma)$vectors %*% (diag(sqrt(eigen(sigma)$values)) * stdev)
  points = t(colMeans(points) + A%*%x)
  points <- data.frame(points)
  colnames(points) <- c("x","y")
  points
}

generate_correlated_data <- function(num_obs, corr){
  desired_correlation <- corr
  correlation <- -20 #intialize the correlation to something you'd never get. 
  
  while(correlation != desired_correlation){
    mu <- rep(0,2)
    Sigma <- matrix(corr, nrow=2, ncol=2) + diag(2)*(1-corr)
    rawvars <- mvrnorm(n=num_obs, mu=mu, Sigma=Sigma)
    
    correlation <- round(cor(rawvars)[2],3)
  }
  rawvars
}


#Create a function akin to the GGI but can be used for high dimensional data. 
PCA_MAD <- function(Subject_Vector, RefPop_Matrix, Min_VAF){
  #Compute the mean and std dev of each variable to use later to center and scale the subject data
  refpop_means <- colMeans(RefPop_Matrix)
  refpop_sd <- apply(RefPop_Matrix,2,sd)
  
  #run the principal component analysis on the RefPop_Matrix.
  p <- prcomp(refpop_sc, center=T, scale=T, retx=T)
  
  #Determine number of PC's to keep based on Min_VAF
  EigenSum <- sum(p$sdev^2)
  EigenValues <- p$sdev^2
  PC_PercentVarExplained <- t(matrix(abs(EigenValues/EigenSum)))
  numPCs<- 1
  repeat{
    if(sum(PC_PercentVarExplained[1:numPCs])>=Min_VAF){
      break
    }
    numPCs<- numPCs+ 1
  }
  
  #Calculate projections
  refpop_projs<- refpop_sc %*% p$rotation[,1:numPCs]
  
  # Get the means of the PCs for the reference group
  # These should all be near zero.
  PCmeans<- t(as.matrix(colMeans(refpop_projs)))
  
  # Get the SDs of the Time Distance PCs. Should = sqrt(Eigenvalues)
  PCsds<- t(as.matrix(apply(refpop_projs,2,sd)))  
  
  Subject_Vector <- data.frame(t(Subject_Vector))
  Subject_Vector <- data.frame(scale(Subject_Vector, center = refpop_means, scale = refpop_sd))
  
  #Convert to a matrix for PC projection
  glpopm_sc<- data.matrix(Subject_Vector)
  
  #Multiply the patient data by the eigenvectors calculated on the normal population. Result is the PC values for the given subject.
  glpopm_projs<- glpopm_sc%*% p$rotation[,1:numPCs]
  
  dist <- (glpopm_projs-PCmeans[,1:numPCs])/PCsds[,1:numPCs]
  abs_dist <- abs(dist)
  sum(abs_dist)/numPCs
}

#https://stackoverflow.com/questions/37773469/r-random-distribution-with-predefined-min-max-mean-and-sd-values
rgbeta <- function(n, mean, var, min = 0, max = 1)
{
  dmin <- mean - min
  dmax <- max - mean
  
  if (dmin <= 0 || dmax <= 0)
  {
    stop(paste("mean must be between min =", min, "and max =", max)) 
  }
  
  if (var >= dmin * dmax)
  {
    stop(paste("var must be less than (mean - min) * (max - mean) =", dmin * dmax))
  }
  
  # mean and variance of the standard beta distributed variable
  mx <- (mean - min) / (max - min)
  vx <- var / (max - min)^2
  
  # find the corresponding alpha-beta parameterization
  a <- ((1 - mx) / vx - 1 / mx) * mx^2
  b <- a * (1 / mx - 1)
  
  # generate standard beta observations and transform
  x <- rbeta(n, a, b)
  y <- (max - min) * x + min
  
  return(y)
}

F_G_Test <- function(input_matrix){
  x <- input_matrix
  nvar <- ncol(x)
  n <- nrow(x)
  R <- cor(x)
  log_e_det <- determinant(R)$modulus[1] 
  Fchi <- -(nrow(x) - 1 - (1/6) * (2 * nvar + 5)) * log_e_det
  df <- 1/2 * (nvar) * (nvar -  1)
  return_list <- list()
  return_list$Test_Stat <- Fchi
  return_list$df <- df
  return_list$p_val <- 1 - pchisq(Fchi, df)
  return_list
}

Dist <- function(Subj, Ref,Min_VAF = 1){
  n_vars <- length(Subj)
  n_obs <- nrow(Ref)
  
  
  ### Non_PCA Distance measures ####################
  refpop_means <- colMeans(Ref)
  refpop_sd <- apply(Ref,2,sd)
  dist <- (Subj-refpop_means)/refpop_sd
  
  #Euclidean
  k=2
  dist_squared <- (abs(dist))^k
  Euclidean_L2 <- sum(dist_squared)^(1/k)
  
  #MAD
  dist_abs <- abs(dist)
  MAD <- sum(dist_abs)/n_vars
  ####################################################
  
  ### PCA Distance measures ##########################
  
  #Compute the mean and std dev of each variable to use later to center and scale the subject data
  # refpop_means <- colMeans(Ref)
  # refpop_sd <- apply(Ref,2,sd)
  refpop_sc <- scale(Ref,scale=T,center=T)
  refpop_cov<- cov(refpop_sc)
  
  #run the principal component analysis on the Ref data
  p <- prcomp(Ref, center=T, scale=T, retx=T)
  
  #Determine number of PC's to keep based on Min_VAF
  EigenSum <- sum(p$sdev^2)
  EigenValues <- p$sdev^2
  PC_VAF <- t(matrix(abs(EigenValues/EigenSum)))
  cum_PC_VAF <- round(cumsum(PC_VAF),5) #rounding is necessary for logical comparison to work.   #https://stackoverflow.com/questions/2769510/numeric-comparison-difficulty-in-r
  numPCs <- min(which(cum_PC_VAF >= Min_VAF))
  
  #Calculate projections
  refpop_projs<- p$x
  
  # Get the means of the PCs for the reference groups of normal subjects
  # These should all be near zero.
  PCmeans<- colMeans(refpop_projs)
  
  # Get the SDs of the Ref PCs. Should = sqrt(Eigenvalues)
  PCsds <- p$sdev
  
  #scale and center Subj based on the mean and sd of Ref. 
  Subj_sc <- (Subj-refpop_means)/refpop_sd
  
  #Multiply the subject data by the eigenvectors calculated on the normal population. Result is the PC values for the given subject.
  Subj_projs<- Subj_sc%*% p$rotation[,1:numPCs]
  
  dist <- (Subj_projs-PCmeans[1:numPCs])/PCsds[1:numPCs]
  
  #Euclidean
  k=2
  dist_squared <- (abs(dist))^k
  Euclidean_L2_PCA <- sum(dist_squared)^(1/k)
  
  #MAD
  dist_abs <- abs(dist)
  MAD_PCA <- sum(dist_abs)/numPCs
  
  ########################################################
  Last_EigVal_GT_1 <- max(which(EigenValues > 1))
  
  
  return_list <- list()
  return_list$Euclidean <- round(Euclidean_L2,5)
  return_list$Euclidean_PCA <- round(Euclidean_L2_PCA,5)
  return_list$MAD <- round(MAD,5)
  return_list$MAD_PCA <- round(MAD_PCA,5)
  return_list$VAF <- cum_PC_VAF[numPCs]
  return_list$num_PCs <- numPCs
  return_list$Last_EigVal_GT_1 <- Last_EigVal_GT_1
  return(return_list)
}

generate_hd_correlated_data <- function(n, p, corr,constant_cov_matrix = T,mean=0){
  mu <- rep(mean,p)
  if(constant_cov_matrix==TRUE){
    Sigma <- matrix(corr, nrow=p, ncol=p) + diag(p)*(1-corr)
  }else{
    variance <- (corr - .001) * (.999 - corr)*.15 #var must be less than (mean - min) * (max - mean)
    Matrix_Values <- rgbeta(p*p, mean = corr, var = variance, min = .001, max = .999) #create matrix with mean = corr, but min=0 and max=1
    # if(all_positive_corr == FALSE){
    #   Matrix_Values <- Matrix_Values*sample(c(1,-1),n*n,replace=T)
    # }
    S <- matrix(Matrix_Values, p, p)
    S <- forceSymmetric(S)
    for(i in 1:p){S[i,i] <- 1} #ones along the diagonal
    Sigma <- S
  }
  rawvars <- mvrnorm(n, mu=mu, Sigma=Sigma,tol=1)
  rawvars
}

simulate_distance <- function(n,p,corr, Min_VAF = 1,Subj_Type = 'All_1s',constant_cov_matrix = T){
  # Return_List <- list() #initialize a list
  Ref <- generate_hd_correlated_data(n,p,corr,constant_cov_matrix) #generate random ref population with standardized variables. 
  if(Subj_Type == 'Random'){Subj <- generate_hd_correlated_data(1,p,corr,constant_cov_matrix)}
  if(Subj_Type == 'Random_1s'){#generate random vector of length p of half +1s and half -1s
    Subj <- sample(c(rep(1,ceiling(p/2)),c(rep(-1,ceiling(p/2)))),p,replace = F) #take the ceiling of p to account for p being odd. 
    } 
  if(Subj_Type == 'All_1s'){Subj <- rep(1,p)}
  
  # Subj <- if(Random_Subj) sample(c(1,-1),n,replace = T) else rep(1,n)
  
  distance <- Dist(Subj,Ref,Min_VAF)
  F_G_Test_Results <- F_G_Test(Ref)
  abs_mean_correlation <- mean(abs(cor(Ref)[lower.tri(cor(Ref), diag = FALSE)]))
  
  
  results <- data.frame(n,p,corr, constant_cov_matrix,Min_VAF,abs_mean_correlation,Subj_Type,VAF = distance$VAF,
                        num_PCs = distance$num_PCs,
                        Euclidean = distance$Euclidean,
                        Euclidean_PCA = distance$Euclidean_PCA,
                        MAD = distance$MAD,
                        MAD_PCA = distance$MAD_PCA,
                        FG_Chi_Sq = F_G_Test_Results$Test_Stat,
                        FG_df = F_G_Test_Results$df,
                        FG_p_val = F_G_Test_Results$p_val,
                        Last_EigVal_GT_1 = distance$Last_EigVal_GT_1
  )
  return(results)
}

repeat_distance_simulation <- function(times, n,p,corr, Min_VAF = 1, Subj_Type = 'All_1s',constant_cov_matrix = T){
  sim_results <-  mapply(simulate_distance,rep(n,times),MoreArgs = list(p,corr, Min_VAF,Subj_Type, constant_cov_matrix),SIMPLIFY=FALSE)
  return(rbindlist(sim_results))
}


