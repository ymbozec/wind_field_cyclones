## wind_field_cyclones
Model to estimate maximum mean wind speed at any distance of a cyclone track.

The script can be used to predict mean wind speed at multiple locations given their distance to past cyclone tracks. Requires data on historical cyclone with name, date, spatial coordinates, central pressure, mean radius of maximum mean wind speed and maximum mean wind speed in the vincinity of the cyclone centre. For every track, cyclone characteristics are interpolated at hourly intervals, and the distance to each location is calculated. The forward motion of the cyclone is estimated between two consecutive positions, as well as the clockwise angle between the cyclone track (for asymetry calculations) and the line separating each reef to the moving eye position. Maximum mean wind speed is then estimated based on the distance between each location and the cyclone centre following Holland (1980)'s wind field model with asymmetry (based on McConochie et al. 2004). For each location, the highest value of mean wind speed is recorded. Although a location may be affected by multiple cyclones within a single season, only the highest value of maximum wind speed across all cyclone events for that season is retained.

The method is used to reconstruct reef exposure to recent cyclones on the Great Barrier Reef (GBR) for the ecological model ReefMod-GBR (v.7.2 and later versions). An extra step for ReefMod is to convert wind speed into a cyclone category defined on the Saffir-Simpson scale (Bozec et al. 2022). Note that the Australian Bureau of Meteorology reports mean wind speed over a 10-minute period, while the Saffir-Simpson scale is based on maximum sustained 1-minute winds. So mean wind speeds need to be converted into 1-min sustained winds equivalent following the relationship: MSW1 = MSW10*1/0.88 (MWO 2017, Kruk et al. 2010).

Requires two datasets: 
1) the spatial coordinates of the locations of interest
2) a data table listing the position and characteristics of cyclones (typically extracted from international databases of best cyclone tracks).

An example is given with the spatial coordinates of 3,806 reef centroids of the GBR (https://github.com/ymbozec/REEFMOD.7.2_GBR) and Australian cyclone tracks for the period 2008-2026 derived from the Australian Tropical Cyclone Database (https://www.bom.gov.au/cyclone/tropical-cyclone-knowledge-centre/databases/).

## Contact
Yves-Marie Bozec, The University of Queensland (y.bozec@uq.edu.au)
