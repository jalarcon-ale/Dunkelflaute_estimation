## Spatial Dependence
To estimate the spatial extent of dunkelflautes, we consider one spatial dependence measure from the literature: the F-madogram(see [3, 8]).
The F-madogram for locations r1 and r2 separated by h = ∥r1 − r2∥ is
defined as:
νF (h) =
1
2
E[|F{Z(r1)} − F{Z(r2)}|] , (3)
where Z(r) represents yearly minima of CFtot at location r, and F is its cumulative
distribution function. Notably, νF (h) = 1/6 for independent year minima
and 0 ≤ νF (h) < 1/6 when dependent [3].


The Spatial_correlation.R file contains a script to estimate the spatial correlation of low energy values corresponding to the energy mix
 

## Dunkelflaute estimation
The Estimate_droughts.R file contains a script to estimate return periods of widespread energy droughts across Europe. This is done with respect to more than fourty years of daily capacity factor data on the ERA5 grid. The PhaseRando.R file contains the functions needed to generate spatial extremes. The method is laid out by Hans Van de Vyver, in "Fast generation of high-dimensional spatial extremes". 
