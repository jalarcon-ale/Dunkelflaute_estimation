library(ncdf4)
library(evd)
library(readr)	#library to read in csv files
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)

###Calling the functions in PhaseRando.R
source("/home/alejala/fileserver/home/script/PhaseRando.R")

###Retriving daily SPV capacity factors
setwd("/mnt/HDS_CORDEXBE_RMIB/hvijver/ReNew/SPV/daily")
spv_nc_file = nc_open("cds-renewables-daily-europe-solar_photovoltaic_power_generation-1979_2018.nc")
###Retriving daily WON capacity factors
setwd("/mnt/HDS_CORDEXBE_RMIB/hvijver/ReNew/WON/daily")
won_nc_file = nc_open("cds-renewables-daily-europe-wind_power_generation_onshore_1979_2018.nc")


###Retrieving data from nc_file
spv_df = ncvar_get(spv_nc_file, varid="spv_cf")
won_df = ncvar_get(won_nc_file, varid="woncfr")
lat = ncvar_get(spv_nc_file, "latitude")
lon = ncvar_get(spv_nc_file, "longitude")
times= ncvar_get(spv_nc_file, "time")
time_steps = length(times)
nc_close(spv_nc_file)
nc_close(won_nc_file)


### Regional data
### Region changes with delta
### DELTA=0.7,...,1
lon_min = 4.25  # Western boundary
lon_max = 9.75   # Eastern boundary
lat_min = 47.25   # Southern boundary
lat_max = 51.75   # Northern boundary

### DELTA=0.6
lon_min = 4.5  # Western boundary
lon_max = 9.50   # Eastern boundary
lat_min = 48   # Southern boundary
lat_max = 51   # Northern boundary

### DELTA=0.5
lon_min = 4.75  # Western boundary
lon_max = 9.25   # Eastern boundary
lat_min = 48   # Southern boundary
lat_max = 51   # Northern boundary

### DELTA=0.4
lon_min = 5  # Western boundary
lon_max = 9   # Eastern boundary
lat_min = 48.25   # Southern boundary
lat_max = 51   # Northern boundary

### DELTA=0.3
lon_min = 5.25  # Western boundary
lon_max = 9   # Eastern boundary
lat_min = 48.25   # Southern boundary
lat_max = 50.75   # Northern boundary

### DELTA=0.2
lon_min = 5.75  # Western boundary
lon_max = 8.25   # Eastern boundary
lat_min = 48.75   # Southern boundary
lat_max = 50.25   # Northern boundary

### DELTA=0.1
lon_min = 6  # Western boundary
lon_max = 8   # Eastern boundary
lat_min = 48.75   # Southern boundary
lat_max = 50.25   # Northern boundary

### DELTA=0
lon_min = 6.25  # Western boundary
lon_max = 7.75   # Eastern boundary
lat_min = 49   # Southern boundary
lat_max = 50   # Northern boundary


### Subsetting spatiotemporal data according to the extent
lon_idx = which(lon >= lon_min & lon <= lon_max)
lat_idx = which(lat >= lat_min & lat <= lat_max)
### Data
West_spv_df = spv_df[lon_idx, lat_idx, ]   
West_won_df = won_df[lon_idx, lat_idx, ]   
lon_West = lon[lon_idx]
lat_West = lat[lat_idx]
num_locations=length(lon_West)*length(lat_West) 


### Data coordinates in 2D
coord.grid=expand.grid(longitude=lon_West, latitude=lat_West)
n.grid=nrow(coord.grid)

###Matrix with distance between gridpoints
D = as.matrix(dist(coord.grid))
D=111*D


### Energy mix
DELTA=0
data_df = (1-DELTA)*West_won_df+DELTA*West_spv_df


### The spatial data is transformed to 2D, where the order in the 
### columns is the one induced by expand.grid
Y=matrix(-data_df, nrow = time_steps, ncol = num_locations,
byrow=T, dimnames=list(time_evolution,locations))

#colnames(Y)


###Transforming data to Pareto form, while keeping its form
P=Pareto_margins(Y)


###To be done once: retrieving only the extreme fields
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

### Saving and loading data
write.csv(estimates,"Estimates_Delta_00.csv", row.names = FALSE)
 


#############################################################

setwd("/home/alejala/fileserver/home/script/plots")


# Generate filenames and read them into a list
estimates_list <- lapply(0:10, function(i) {
  delta_str <- if (i == 10) "1" else sprintf("0%d", i)
  read_csv(paste0("Estimates_Delta_", delta_str, ".csv"))
})



# Create a data frame for plotting
plot_data <- data.frame(
  Delta = seq(0, 1, by = 0.1),
  ### Extract the return period for 100% joint exceedance
  seventy = sapply(estimates_list, function(df) df[[3, 3]]),
  eighty = sapply(estimates_list, function(df) df[[4, 3]]),
  ninety = sapply(estimates_list, function(df) df[[5, 3]]),
  all = sapply(estimates_list, function(df) df[[6, 3]]) 
)




# Convert to long format
plot_data_long <- pivot_longer(
  plot_data,
  cols = -Delta,
  names_to = "Joint_exceedance",
  values_to = "Return_Period"
)

# Rename exceedance levels
plot_data_long <- plot_data_long %>%
  mutate(Joint_exceedance = recode(Joint_exceedance,
    seventy = "70%",
    eighty = "80%",
    ninety = "90%",
    all = "100%"
  ))

# Ensure correct order
plot_data_long$Joint_exceedance <- factor(
  plot_data_long$Joint_exceedance,
  levels = c("70%", "80%", "90%", "100%")
)

# Plot
p <- ggplot(plot_data_long, aes(x = Delta, y = Return_Period, fill = Joint_exceedance)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_log10() +
  scale_fill_brewer(palette = "Set2") +  # Colorblind-friendly and distinguishable
  theme_minimal() +
  labs(
    x = "Delta",
    y = "Return period in years (log scale)",
    fill = "Joint exceedance"
  )

print(p)


  
# Save to PNG
ggsave("estimates_plot.png", plot = p, width = 8, height = 6, dpi = 300)
