This respository supports the project $\textbf{Statistical Modeling of Widespread Energy Droughts}$. The written report can here be found under the name technische_rapport.pdf. Two main calculations are addressed in this project: 
* The estimation of the spatial extent of energy droughts.
* The estimation of return periods of widespread energy droughts.

These estimations can be done with respect to more than fourty years of daily capacity factor data on the ERA5 grid, as available from the Copernicus climate change service: https://cds.climate.copernicus.eu/datasets/sis-energy-derived-reanalysis?tab=overview 

## Spatial Dependence
To estimate the spatial extent of energy droughts, we consider one spatial dependence measure from the literature: the F-madogram (see A. Gobin and H. Van de Vyver, 2021).
The $F$-madogram for locations $r_1$ and $r_2$ separated by $h = ∥r_1 − r_2∥$ is defined as:

$$ v_{F} (h) = \dfrac{1}{2} E[|F{Z(r1)} − F{Z(r2)}|], $$

where $Z(r)$ represents yearly minima of $CF_{tot}$ at location $r$, and $F$ is its cumulative distribution function. Notably, $ν_F (h) = 1/6$ for independent year minima, and $0 ≤ ν_F (h) < 1/6$ when dependent [3].


The Spatial_correlation.R file contains a script to estimate the spatial dependence of extremely low energy values $CF_{tot}$.
 

## Dunkelflaute estimation
The Estimate_droughts.R file contains a script to estimate return periods of widespread energy droughts. The PhaseRando.R file contains the functions needed to generate spatial extremes. The method is laid out by H. Van de Vyver, in "Fast generation of high-dimensional spatial extremes". 
