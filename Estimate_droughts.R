library(ncdf4)
library(evd)
library(readr)	#library to read in csv files
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)

install.packages(c("ncdf4","evd","readr", "ggplot2",  "dplyr", "RColorBrewer" , "tidyr"))


###Calling the functions in PhaseRando.R
source("path/PhaseRando.R")
### Loading in
setwd("path")
### Daily SPV capacity factors
spv_nc_file = nc_open("cds-renewables-daily-europe-solar_photovoltaic_power_generation-1979_2018.nc")
### Daily WON capacity factors
won_nc_file = nc_open("cds-renewables-daily-europe-wind_power_generation_onshore_1979_2018.nc")

###Retrieving data from nc_file
spv_df = ncvar_get(spv_nc_file, varid="spv_cf")
won_df = ncvar_get(won_nc_file, varid="woncfr")
lat = ncvar_get(spv_nc_file, "latitude")
lon = ncvar_get(spv_nc_file, "longitude")
times= ncvar_get(spv_nc_file, "time")
time_steps = length(times)
###By closing the nc files, memory frees up
nc_close(spv_nc_file)
nc_close(won_nc_file)

### Initializing a list of data frames
de lijst = list()

### Region changes with delta
for (DELTA in seq(from = 0, to = 1, by = 0.1)) {
	if (DELTA==0) {
		### DELTA=0
        lon_min = 6.25  # Western boundary
        lon_max = 7.75   # Eastern boundary
        lat_min = 49   # Southern boundary
        lat_max = 50   # Northern boundary
    } 
	
	else if (DELTA==0.1) {
	    lon_min = 6  # Western boundary
        lon_max = 8   # Eastern boundary
        lat_min = 48.75   # Southern boundary
        lat_max = 50.25   # Northern boundary
    }
	
	else if (DELTA==0.2) {
	    lon_min = 5.75  # Western boundary
        lon_max = 8.25   # Eastern boundary
        lat_min = 48.75   # Southern boundary
        lat_max = 50.25   # Northern boundary
	}
	
	else if (DELTA==0.3) {
	    lon_min = 5.25  # Western boundary
        lon_max = 9   # Eastern boundary
        lat_min = 48.25   # Southern boundary
        lat_max = 50.75   # Northern boundary	
	}
	
	else if (DELTA==0.4) {
		lon_min = 5  # Western boundary
        lon_max = 9   # Eastern boundary
        lat_min = 48.25   # Southern boundary
        lat_max = 51   # Northern boundary
    }
	
	else if (DELTA==0.5) {
	    lon_min = 4.75  # Western boundary
        lon_max = 9.25   # Eastern boundary
        lat_min = 48   # Southern boundary
        lat_max = 51   # Northern boundary
	}
	
	else if (DELTA==0.6) {
	    lon_min = 4.5  # Western boundary
        lon_max = 9.50   # Eastern boundary
        lat_min = 48   # Southern boundary
        lat_max = 51   # Northern boundary
    }
	
	else if (DELTA>=0.7) {
		lon_min = 4.25  # Western boundary
        lon_max = 9.75   # Eastern boundary
        lat_min = 47.25   # Southern boundary
        lat_max = 51.75   # Northern boundary
    }

    ###Subsetting spatiotemporal data according to rectangular extent
    lon_idx = which(lon >= lon_min & lon <= lon_max)
    lat_idx = which(lat >= lat_min & lat <= lat_max)
    West_spv_df = spv_df[lon_idx, lat_idx, ]   
    West_won_df = won_df[lon_idx, lat_idx, ]   
    
    ###Retriving the coordinates
    lon_West = lon[lon_idx]
    lat_West = lat[lat_idx]
    num_locations=length(lon_West)*length(lat_West) 
    ### Data coordinates in 2D
    coord.grid=expand.grid(longitude=lon_West, latitude=lat_West)
    n.grid=nrow(coord.grid)
    ###Matrix with distance between gridpoints
    D = as.matrix(dist(coord.grid))
    D=111*D ###roughly changing to km

    ### Energy mix
    data_df = (1-DELTA)*West_won_df+DELTA*West_spv_df

    ### The spatial data is transformed to 2D, where the order in the 
    ### columns is the one induced by expand.grid
    Y=matrix(-data_df, nrow = time_steps, ncol = num_locations, byrow=T)

    ###Transforming data to Pareto form
    P=Pareto_margins(Y)

    ###To be done once: 
    ###Filtering extreme fields
    M_list=EXTREMEData(P)
    M=M_list[[1]]   		###These are the extreme events
    r_O=M_list[[2]]
    days_idx=M_list[[3]]   	### Days with extremely low CF values
    n.exc <- dim(M)[1]   	### Number of extreme fields

    ###Monte Carlo simulations
    n.sim=11000
    ###Creating list of phase-randomised spatial extremes
    sim_list=list(M)
    ###Setting seed for reproductibility
    semilla=1
    for (i in seq(n.sim)) {
	    print(i)
	    sim_list=append(sim_list,list(PRSim(M,r_O,ALPHA,semilla)))
	    semilla=semilla+1
    }

    n.fields=n.exc*(n.sim+1)   ###Number of fields created
    
    ###Using the theory to find tau-quantile
    v_O=r_O/n.grid   ###Lower bound for concurrent extremes

    j=1
    estimates=replicate(5,0)
    ###Looping over proportions 
    for (i in c(0.5,0.6,0.7,0.8,0.9,1)){
	    print(i)
	    estimates[j]=concurrent_optimized(sim_list,n.sim,n.exc,v_O,i)
        j=j+1
    }
    tau_O=1-v_O^(-ALPHA)   ###Using the theory to find tau-quantile
    Joint_exceedance=estimates*(1-tau_O)
    Joint_exceedance=(0.01)*(v_O^ALPHA)*Joint_exceedance

    ###Final data frame 
    estimates=data.frame(
             Proportion_of_locations=c("50%","60%","70%","80%","90%","100%"),
             Joint_exceedance,
             Return_period=1/(365*Joint_exceedance)
             )
    
    ###Appending to the list of data frames
    de_lijst=append(de_lijst,estimates)
	
}





