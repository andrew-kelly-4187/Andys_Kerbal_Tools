# Solving Lambert's problem

# The first attempt will just be Kerbin to Duna as a test
source("KSPParameters.R")
source("Astrodynamics functions.R")

# Let's just pick a certain time in the future
time = 10000000

# Determine the true anomaly for both Kerbin and Duna at this new time
MeanAnomaly_Kerbin = kerbin_meananomaly + (2*pi*time)/(kerbin_siderealperiod)
MeanAnomaly_Duna = duna_meananomaly + (2*pi*time)/(duna_siderealperiod)

# Now calculate eccentric anomaly for both
EccentricAnomaly_Kerbin = fn_kepler(MeanAnomaly_Kerbin, kerbin_eccentricity)
EccentricAnomaly_Duna = fn_kepler(MeanAnomaly_Duna, duna_eccentricity)

# Now calculate true anomaly for both
TrueAnomaly_Kerbin = 2 * atan(tan(0.5 * EccentricAnomaly_Kerbin)*sqrt((1 + kerbin_eccentricity)/(1 - kerbin_eccentricity)))
TrueAnomaly_Duna = 2 * atan(tan(0.5 * EccentricAnomaly_Duna)*sqrt((1 + duna_eccentricity)/(1 - duna_eccentricity)))

# Now convert to perifocal coordinates
r_Kerbin = (kerbin_angularmomentum^2)/(kerbol_mu*(1+kerbin_eccentricity*cos(TrueAnomaly_Kerbin)))*c(cos(TrueAnomaly_Kerbin),sin(TrueAnomaly_Kerbin),0)
v_Kerbin = (kerbol_mu/kerbin_angularmomentum)*c(-sin(TrueAnomaly_Kerbin), (kerbin_eccentricity + TrueAnomaly_Kerbin), 0)

r_Duna = (duna_angularmomentum^2)/(kerbol_mu*(1+duna_eccentricity*cos(TrueAnomaly_Duna)))*c(cos(TrueAnomaly_Duna),sin(TrueAnomaly_Duna),0)
v_Duna = (kerbol_mu/duna_angularmomentum)*c(-sin(TrueAnomaly_Duna), (duna_eccentricity + TrueAnomaly_Duna), 0)

