


##########################################
###Function to retrieve the extreme fields
###Input is the spatiotemporal dataset 
###transformed to Pareto margins

EXTData<-function(Z) {
	###Z is a spatio-temporal data set with Pareto margins
    time_steps=dim(Z)[1]
    num_locations=dim(Z)[2]
    ###The spatial sum is the radial component
    spatial_sum=rep(0,time_steps)  ### sum per day
    for (i in 1:time_steps){
	    spatial_sum[i]=sum(Z[i,])
    }
    ###Indixes corresponding to extreme days  
    r_O=quantile(spatial_sum,probs=0.95)
    days_idx = which(spatial_sum>r_O)
    ext_days = length(days_idx)   ###Number of extreme fields

    ###Retrieving extreme events
    X=array(0,dim=c(ext_days,num_locations))
    k=1
    for (i in days_idx) {	
        X[k,]=Z[i,]
        k=k+1
    }
    out=list(X,r_O)
    return (out)
}




#######################################################
#############FUNCTION for phase randomisation 
####(input is matrix M; rows: time, columns:location)


PRSim<-function(M,r_O,ALPHA) {
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


