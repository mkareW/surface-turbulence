# surface-turbulence
R scripts for calculation of surface turbulent fluxes over snow/ice

mobulk2: older version 

mobulk3: R function to calculate surface turbulent fluxes (latent and sensible heat fluxes), with the "bulk aerodynamic" method. You can use either Brutsaert or Monin-obukhov stability functions.
Input data from a meteorological station located above the ice.

Calculates the error on the fluxes due to uncertainties on the measurements.

The function takes as entry:

#V windspeed in ms-1

#Ta air temperature in K

#RH air relative humidity in %

#p air barometric pressure in mbars

#Zv,Zt,Zrh : height of the measurements (V,T,RH) in m 

#Ts: temperature of the surface in K

#z0,z0t,z0q: surface roughness lengths in m
