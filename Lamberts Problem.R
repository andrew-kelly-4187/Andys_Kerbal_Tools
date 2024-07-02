# Solving Lambert's problem

# The first attempt will just be Kerbin to Duna as a test
source("KSPParameters.R")
source("Astrodynamics functions.R")

# Let's just pick a certain time in the future
time = 1000000

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

# Let's use the method in Prussing and Conway "Orbital Mechanics"
# We will also work in canonical units
c_vec = (r_Duna - r_Kerbin)/canonical_DU
c = fn_veclength(c_vec)
r_1 = fn_veclength(r_Kerbin)/canonical_DU
r_2 = fn_veclength(r_Duna)/canonical_DU
s = 0.5*(r_1 + r_2 + c)

phase = acos((r_1^2+r_2^2-c^2)/(2*r_1*r_2))

a_m = s/2 # Minimum semimajor axis of transfer, DU
beta_m = 2 * asin(sqrt((s-c)/s))
t_m = sqrt((s^3)/(8))*(pi - beta_m + sin(beta_m)) # Transfer time for minimum semimajor axis orbit, TU


# Now we bound the transfer time:
t_p = (sqrt(2)/3)*(s^(3/2) - sign(sin(phase)) *(s-c)^(3/2))
# This is the minimum possible transfer time

# Now we provide the transfer time we want to determine the semimajor axis for:
tof = 10 #TU

# Now we solve for a using an interval bisection approach
f <- function(a,s,c,t){
  alpha = 2*asin(sqrt(s/(2*a)))
  beta = 2*asin(sqrt((s-c)/(2*a)))
  value = alpha^(3/2) *(alpha - beta - (sin(alpha) - sin(beta))) - t
  return(value)
}

# Initial interval
a1 = a_m
a2 = 1.5*a_m
step = a2-a1
while(step > 0.001){
  
  a_mid = 0.5*(a1 + a2)
  val_mid = f(a_mid,s,c,tof)
  
  if(val_mid <0){
    a1 = a1
    a2 = a_mid
  } else {
    a1 = a_mid
    a2 = a2
  }

  if(a2-a1 > step){stop('Algorithm is diverging')}
  step = a2-a1

  a_root = 0.5*(a1 + a2)
}

# Work out alpha and beta, considering the quadrants and branches
alpha = 2*asin(sqrt(s/(2*a)))
beta = 2*asin(sqrt((s-c)/(2*a)))
if(tof > t_m){alpha = 2*pi - alpha}
if(phase >= pi){beta = -beta}

# Now to determine the terminal velocity vectors
A = pracma::cot(alpha/2)*sqrt(1/(4*a))
B = pracma::cot(beta/2)*sqrt(1/(4*a))
u1 = r_Kerbin/fn_veclength(r_Kerbin)
u2 = r_Duna/fn_veclength(r_Duna)
uc = c_vec/c

v1 = (B+A)*uc - (B-A)*u1
v2 = (B+A)*uc - (B-A)*u2

# Delta-V changes required
dv1 = abs(fn_veclength(v1 - (v_Kerbin)/canonical_VU))
dv2 = abs(fn_veclength(v2 - (v_Duna)/canonical_VU))
dv_total <- (dv1 + dv2) * canonical_VU
