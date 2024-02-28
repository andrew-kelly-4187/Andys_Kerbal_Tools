# In this file I will be exploring Python
# I always find I learn better with a project, so I will be building a series
# of tools and functions for my hobby/passion: astrodynamics, orbital mechanics
# and mission planning.

#### Imports
import numpy as np
import math as math

# First is a function for taking the cartesian position and velocity vectors
# of a body, and returning the orbital parameters.

def fn_cartesian2elements(r_vec,v_vec, mu_val):
    '''Return the Keplerian elements and orbital parameters from position and velocity data.

    Parameters
    -----------
    r_vec: A three element vector for the position of the body.
    v_vec: A three element vector for the velocity of the body.
    mu_val: The gravitational parameter for the centrally orbited body

    Returns:
    -----------
    Vectors:
        h_vec: The angular momentum vector for the orbit.
        e_vec: The eccentricity vector for the orbit.
        n_vec: The vector of nodes for the orbit.
    Scalars:
        h_val: The angular momentum for the orbit.
        eccentricity_val: The eccentricity for the orbit.
        energy_val: The specific mechanical energy for the orbit.
        a_val: The semimajor axis for the orbit.
        parameter_val: The parameter (semilatus rectum of the orbit).
        inclination_val: The inclination of the orbit.
        LAN_val: The Longitude of the Ascending Node of the orbit.
        ArgOfPeri_val: The argument of periapsis for the orbit.
        TrueAnomaly_val: The true anomaly along the orbital path.
        FlightAngle_val: The angle between the normal plane and prograde.
        orbit_character: A string describing the orbit conic section.
    '''

    # First check that the entries are vectors
    if isinstance(r_vec, (list, tuple, np.ndarray)) == False:
        print('Error: position vector r_vec must be an array of length three')
        return
    elif len(r_vec) != 3:
        print('Error: position vector r_vec must be an array of length three')
        return
    
    if isinstance(v_vec, (list, tuple, np.ndarray)) == False:
        print('Error: position vector r_vec must be an array of length three')
        return
    elif len(v_vec) != 3:
        print('Error: position vector r_vec must be an array of length three')
        return
    
    # Convert vector inputs to numpy arrays
    r_vec = np.array(r_vec)
    v_vec = np.array(v_vec)

    # Get scalar r, v
    r_val = np.linalg.norm(r_vec)
    v_val = np.linalg.norm(v_vec)

    # Return h vector
    h_vec = np.cross(r_vec,v_vec)
    h_val = np.linalg.norm(h_vec)

    # Return n vector
    n_vec = np.cross([0,0,1],h_vec)
    n_val = np.linalg.norm(n_vec)

    # Return e vector   
    e_vec = ((v_val**2 - mu_val/r_val)*r_vec - np.dot(r_vec,v_vec)*v_vec)/mu_val
    eccentricity_val = np.linalg.norm(e_vec)

    # Return specific mechanical energy
    energy_val = v_val**2/2 - mu_val/r_val

    # Return semimajor axis
    a_val = -mu_val/(2*energy_val)

    # Return semilatus rectum
    p_val = h_val**2/mu_val

    # Return inclination
    inclination_val = math.degrees(math.acos(h_vec[2]/h_val))

    # Return LAN
    LAN_val = math.degrees(math.acos(n_vec[0]/n_val))
    if n_vec[1] > 0 and LAN_val > 180: # Quadrant ambiguity check
        LAN_val = 360 - LAN_val

    # Return argument of periapsis
    ArgOfPeri_val = math.degrees(math.acos(np.dot(n_vec, e_vec)/(n_val * eccentricity_val)))
    if e_vec[2] > 0 and ArgOfPeri_val > 180:
        ArgOfPeri_val = 360 - ArgOfPeri_val

    # Return true anomaly
    TrueAnomaly_val = math.degrees(math.acos(np.dot(e_vec, r_vec)/(r_val * eccentricity_val)))
    if (np.dot(r_vec,v_vec)) > 0 and TrueAnomaly_val > 180:
        TrueAnomaly_val = 360 - TrueAnomaly_val

    # Return flight angle
    FlightAngle_val = math.degrees(math.acos(h_val/(r_val * v_val)))
    if np.sign(np.dot(r_vec,v_vec)) != np.sign(FlightAngle_val):
        FlightAngle_val = 360 - FlightAngle_val

    # Return the character of the orbit
    if 1/a_val > 0:
        if 0 <= eccentricity_val <1:
            orbit_character = "Ellipse"
        elif eccentricity_val == 1:
            orbit_character = "Rectilinear ellipse"
    elif 1/a_val == 0:
        orbit_character = "Parabola"
    elif 1/a_val < 0:
        if eccentricity_val == 1:
            orbit_character = "Rectilinear hyperbola"
        elif eccentricity_val > 1:
            orbit_character = "Hyperbola"

    print("\nOrbit type:\n",orbit_character)
    print("\nSemimajor axis:\n",a_val)

    return {
        'Vector of momentum:':h_vec,
        'Vector of eccentricity:':e_vec, 
        'Vector of nodes:':n_vec,
        'Eccentricity:':eccentricity_val,
        'Specific mechanical energy:':energy_val,
        'Semi-major axis:':a_val,
        'Parameter:':p_val,
        'Inclination:':inclination_val,
        'Longitude of the ascending node:':LAN_val,
        'Argument of periapsis:':ArgOfPeri_val,
        'True anomaly:':TrueAnomaly_val,
        'Flight angle:':FlightAngle_val,
        'Orbital character:':orbit_character}

r = [8752,3228,214] #km
v = [-3.2,1.2,-2.2] #km/s
mu = 3.986e5 #km^3/s^2

a = fn_cartesian2elements(r,v,mu)

# This is the notation to extract from a dictionary.
print(a)
print(a.get('Eccentricity'))

