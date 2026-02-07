This respository supports the project $\textbf{Statistical Modeling of Widespread Energy Droughts}$, which was worked out during an internship at the RMI as part of postgraduate studies. The written report can be found in this repository under the name Technisch_rapport.pdf. Two main calculations are addressed in this project: 
* The estimation of the spatial extent of energy droughts.
* The estimation of return periods of widespread energy droughts.

These estimations can be done with respect to energy indicators for Europe derived from ERA5-reanalysis (1979-present), as available from the Copernicus climate change service: https://cds.climate.copernicus.eu/datasets/sis-energy-derived-reanalysis?tab=overview. There, we selected the capacity factors of (i) Solar photovoltaic power generation (SPV), and (ii) Onshore wind power generation (WON) and stored the NetCDF data (daily resolution) on the RMI-server.

### Spatial Dependence
To estimate the spatial extent of energy droughts, we consider one spatial dependence measure from the literature: the F-madogram (see A. Gobin and H. Van de Vyver, 2021). The Spatial_correlation.R file contains a script to estimate such spatial dependence.
 
### Dunkelflaute estimation
The PhaseRando.R file contains the functions needed to generate synthetic spatial extremes, which are in turn used to estimate return periods of widespread energy droughts, as done in the Estimate_droughts.R file. The method followed is laid out by H. Van de Vyver, 2024.
