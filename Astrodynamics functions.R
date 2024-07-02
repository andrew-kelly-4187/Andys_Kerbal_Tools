#------------------------------------------------------------------#
# List of functions that are used in astrodynamics
#------------------------------------------------------------------#

# This script contains functions used in other scripts, be they basic
# maths functions that aren't provided elsewhere, or specific astro
# functions

# As convention, all are prefixed "fn_"

#------------------------------------------------------------------#
# Vector functions                                              ####
#------------------------------------------------------------------#

#### Cartesian length of a vector
# because length() is the number of elements
fn_veclength <- function(a){

  l = sqrt(a %*% a)
  return(l[1])
  
}

#### Cross Product of two 3D vectors
fn_crossprod <- function(a, b) {
  if(length(a)!=3 || length(b)!=3){
    stop("Cross product is only defined for 3D vectors.");
  }
  i1 <- c(2,3,1)
  i2 <- c(3,1,2)
  return (a[i1]*b[i2] - a[i2]*b[i1])
}

#### Dot Product of two vectors
fn_dotprod <- function(a,b){
  
  s = a %*% b
  return(s[1])
  
}

#------------------------------------------------------------------#
# Numerical method functions                                    ####
#------------------------------------------------------------------#

#### Newton Raphson solver for Kepler's equation
fn_kepler <- function(M, e){
  # Initial guess at E
  E = M
  
  h = (E-abs(e)*sin(E)-M) / (1-abs(e)*cos(E))
  
  i=0 # Counter to escape if the loop doesn't terminate
  while (abs(h) >= 0.0001){
    E = E - h
    h = (E-abs(e)*sin(E)-M) / (1-abs(e)*cos(E))
    i = i+1
    if(i>40){
      cat("ERROR: Solver failed to converge")
      break
    }
  }
  return(E)
}

#------------------------------------------------------------------#
# Astrodynamics functions                                       ####
#------------------------------------------------------------------#

#### findxy
# Computes all x, y for single and multi-revolution solutions




#### Lambert's problem solver
# This is the Izzo algoritm
# https://www.esa.int/gsp/ACT/doc/MAD/pub/ACT-RPR-MAD-2014-RevisitingLambertProblem.pdf

fn_Lambert <- function(r1_vec, r2_vec, t, mu){
  
  # Check for incorrect inputs
  if(length(r1_vec) != 3 | length(r2_vec) != 3){
    stop("ERROR: r_1 and r_2 must be vectors of length 3")
  }
  
  if(t < 0 | mu < 0){
    stop("ERROR: t and mu must be greater than zero")
  }
  
  c_vec = r2_vec - r1_vec
  c = fn_veclength(c_vec)
  r1 = fn_veclength(r1_vec)
  r2 = fn_veclength(r2_vec)
  s = 0.5 * (r1 + r2 + c)
  i_hat_r1 = r1_vec/r1
  i_hat_r2 = r2_vec/r2
  i_hat_h = fn_crossprod(i_hat_r1, i_hat_r2)
  lambda = sqrt(1-c/s)
  
  if(r1_vec[1]*r2_vec[2]-r1_vec[2]*r2_vec[1] <0){
    lambda = -lambda
    i_hat_t1 = fn_crossprod(i_hat_r1, i_hat_h)
    i_hat_t2 = fn_crossprod(i_hat_r2, i_hat_r2)
  } else {
    i_hat_t1 = fn_crossprod(i_hat_h, i_hat_r1)
    i_hat_t2 = fn_crossprod(i_hat_h, i_hat_h)
  }
  TimeT = t * sqrt(2*mu/(s^3))
  
  
}
  


