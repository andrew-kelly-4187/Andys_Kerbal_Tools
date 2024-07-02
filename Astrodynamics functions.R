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
# Astrodynamics functions                                       ####
#------------------------------------------------------------------#

# 

