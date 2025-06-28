##########################################
###Function to retrieve declustered extreme fields
###Input is the spatiotemporal dataset transformed to Pareto margins
###Run_length is the minimum number of non-extreme observations between clusters of exceedances
###The output is a list with three components: 
###the extreme-fields
###the 99th percentile of the spatial sum 
###the indices in the dataset that correspond to the extreme fields

EXTREMEData <- function(Z, run_length = 7) {
  # Z is a spatio-temporal data set with Pareto margins
  time_steps <- dim(Z)[1]
  num_locations <- dim(Z)[2]
  
  # Compute spatial sum (radial component)
  spatial_sum <- rowSums(Z)
  
  # Threshold for extremes (1st percentile = 0.99 quantile)
  r_O <- quantile(spatial_sum, probs = 0.99)
  extreme_idx <- which(spatial_sum > r_O)
  
  # Declustering: Keep only one (maximum) exceedance per run_length window
  declustered_idx <- c()
  i <- 1
  while (i <= length(extreme_idx)) {
    # Start of current cluster
    start_idx <- extreme_idx[i]
    
    # Find indices within run_length days
    cluster_window <- which((extreme_idx - start_idx) <= run_length & (extreme_idx - start_idx) >= 0)
    cluster_idxs <- extreme_idx[cluster_window]
    
    # Keep the index of the maximum spatial sum within the cluster
    max_idx <- cluster_idxs[which.max(spatial_sum[cluster_idxs])]
    declustered_idx <- c(declustered_idx, max_idx)
    
    # Move to first exceedance outside the cluster window
    i <- max(cluster_window) + 1
  }
  
  # Retrieve the extreme fields after declustering
  X <- Z[declustered_idx, , drop = FALSE]  # drop = FALSE preserves matrix structure
	
  out <- list(X = X, r_O = r_O, selected_days = declustered_idx)
  return(out)
}




#######################################################
#############FUNCTION for phase randomisation 
####(input is matrix M of extreme fields; rows: time, columns:location)


PRSim<-function(M,r_O,ALPHA,semilla) {
	n.exc <- dim(M)[1]   ### Number of extreme fields
	n.grid <- dim(M)[2]  ### Number of grid-points
    #################################
    ###I.Perform Fourier Transform
    #################################
    ft<-apply(M,2,fft)          ### Fast Fourier Transfrom	
    modulus <- apply(ft,2,Mod)  ### Extract Moduli
    arg<-apply(ft,2,Arg)        ### Extract Arguments
    #################################
    first_part <- 2:(floor(n.exc/2)+1)
    second_part <- (n.exc+1)-(1:floor(n.exc/2))
    arg_first_part <- arg[first_part, ]
    modulus_first_part <- modulus[first_part, ]
    #################################
    ###II.Consider random phases
    #################################
    set.seed(semilla)
    phases_random <- runif(n=length(first_part), min=-pi, max=pi)
    if (n.exc%%2==0) phases_random[length(first_part)] <-0
    ############################################################
    ###III.Add these random phases to the phases of the original
    ###spatial data
    ############################################################
    ft_prsim <- matrix(NA,n.exc,n.grid)
    ft_prsim[1,] <- ft[1,]
    ft_prsim[first_part, ] <- sapply(1:n.grid, function(i){ complex(
	modulus=modulus_first_part[,i], argument=phases_random+
	arg_first_part[ ,i]) } )
	ft_prsim[second_part,] <- Conj(ft_prsim[first_part, ])
    ############################################################
    ###IV.Perform inverse transform
    ft_prsim_inv <- apply(ft_prsim,2,fft,inverse=TRUE)/n.exc
    M_prsim <- Re(ft_prsim_inv)
    ############################################################
    ###
    r=rep(0,n.exc)
    theta=array(0,dim=c(n.exc,n.grid))    
    for (k in seq(n.exc)) {
	    r[k]=sum(M_prsim[k,])   ###Compute radius
        theta[k,]=M_prsim[k,]/r[k]   ###Compute angular component
        r_star=rgpd(1, loc=r_O, scale=r_O/ALPHA, shape=1/ALPHA)
        M_prsim[k,]=r_star*theta[k,]
    }

    return(M_prsim)
}


######################################################
###Function to estimate the proportion of joint exceedance

###Both functions do the same. The optimized version is much faster
concurrent_optimized <- function(sim_list, n.sim, n.exc, v_O, prop) {
  count_matrix <- sapply(seq_len(n.sim + 1), function(i) {
    rowMeans(sim_list[[i]][1:n.exc, ] > v_O) >= prop
  })
  suma <- sum(count_matrix)
  n.fields <- n.exc * (n.sim + 1)
  estimate <- suma / n.fields
  return(estimate)
}
