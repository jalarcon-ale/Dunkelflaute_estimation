library(ncdf4)
library(SpatialExtremes)
library(latex2exp)

setwd("C:/Users/ajag1/Documents/KMI/Data")
nc_file <- nc_open("cds-renewables-daily-europe-wind_power_generation_onshore_1979_2018_Amin.nc")

#Subsetting from nc_file
won_df = ncvar_get(nc_file, varid="woncfr")
lat = ncvar_get(nc_file, "latitude")
lon = ncvar_get(nc_file, "longitude")
times= ncvar_get(nc_file, "time")
time_steps = length(times)

#Geographical area under study
lon_min = 4  # Western boundary
lon_max = 10   # Eastern boundary
lat_min = 47   # Southern boundary
lat_max = 52   # Northern boundary

lon_idx = which(lon >= lon_min & lon <= lon_max)
lat_idx = which(lat >= lat_min & lat <= lat_max)

West_won_df = won_df[lon_idx, lat_idx, ]  
nc_close(nc_file)
lon_West = lon[lon_idx]
lat_West = lat[lat_idx]
num_locations = length(lon_West)*length(lat_West)

###########################################
#######Changing zero values 

minimum=min(West_won_df[West_won_df != 0 & !is.na(West_won_df)])

#Indexes
zero_idx=which(West_won_df==0)
for (idx in zero_idx){
	West_won_df[idx]=runif(1, min = minimum/10, max = minimum)
}

#Each column should correspond to one location
Y = matrix(-West_won_df, nrow = time_steps, ncol = num_locations,
byrow=T)

#Each row corresponds to one location
coord.grid=expand.grid(longitude=lon_West, latitude=lat_West)
coord.grid=111*coord.grid
n.grid=nrow(coord.grid)


#Retrieving the block minima of SPV for Europe
#setwd("/mnt/HDS_CORDEXBE_RMIB/hvijver/ReNew/SPV/daily")
nc_file <- nc_open("cds-renewables-daily-europe-solar_photovoltaic_power_generation-1979_2018_Amin.nc")
spv_df <- ncvar_get(nc_file, varid="Amin")
West_spv_df = spv_df[lon_idx, lat_idx, ]

##################################################################
###Version with quantile-based bins
###Each bin has (approximately) the same number of points
setwd("/Users/ajag1/Documents/KMI/script/plots")

# Define number of bins 
nbins <- 20      
threshold <- 0.15
distance_at_threshold=replicate(11,0)
k=1

for (DELTA in seq(0,1,0.1)){
    filename = paste0("fmadogram", sep="_", "DELTA_", DELTA, ".png")
	data_df = (1-DELTA)*West_won_df+DELTA*West_spv_df
	#Each column should correspond to one location
	data_matrix = matrix(-data_df, nrow= time_steps, ncol=num_locations,
	byrow=T)
    # Compute the F-madogram
    fmadog = fmadogram(data_matrix,as.matrix(coord.grid),which=c("mado"))
    fmadog_df <- data.frame(dist = fmadog[,1], madogram = fmadog[,2])
    # Create quantile-based breaks
    statD <- hist(fmadog_df$dist, breaks = quantile(fmadog_df$dist, (0:nbins)/nbins), plot = FALSE)
    breaks <- statD$breaks
    # Assign bins using quantile breaks
    fmadog_df$bin <- cut(fmadog_df$dist, breaks = breaks, include.lowest = TRUE)
    # Compute average madogram in each bin
    binned <- aggregate(madogram ~ bin, data = fmadog_df, FUN = mean)
    # Also get bin centers for plotting
    bin_centers <- aggregate(dist ~ bin, data = fmadog_df, FUN = mean)
    # Combine into one data frame
    binned$center <- bin_centers$dist
    # Find the first bin where average madogram reaches or exceeds 0.15
    idx_thresh <- which(binned$madogram >= threshold)[1]
    # Get the corresponding bin center distance
    distance_at_threshold[k] <- binned$center[idx_thresh]
    k=k+1
}

###Creating data frame with colmuns delta and distance at correlation threshold
Dunkelflautes_df=data.frame(Delta=seq(0,1,0.1), Extent=distance_at_threshold)
write.csv(Dunkelflautes_df,"20_bins_Dunkel_Extent.csv", row.names = FALSE)

###Plot of the distance at threshold vs delta value
png("dunkel_spatial_extent_20_bins.png")
plot(seq(0,1,0.1),distance_at_threshold,type="b", col="blue4",
    xlab=TeX('$\\delta$'), ylab = "Distance (km)", main="Spatial extents of Dunkelflautes")
dev.off()
