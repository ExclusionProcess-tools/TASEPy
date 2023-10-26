###############################################################################
#
#       TASEPy v1.0
#                               October 2023
#
#       Based on Ciandrini et al., 2023
#
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from scipy.special import comb as choose

###############################################################################
#
# DEFINE: Returns the maximum number of particles that can fit on the lattice.
#
###############################################################################

def N_max(L, ll=1):
  '''
  Returns the maximum number of particles that can fit onto the lattice of 
  size L. Particle size ll is optional (=1 by default).
  '''
  
  Nmax = 0
  i = 1
  xi = 1
  while xi <= L:
    Nmax += 1
    i += 1
    xi = 1 + (i-1)*ll
    
  return Nmax



###############################################################################
#
# DEFINE: Returns the stacked particle configuration for a given PSA order.
#
###############################################################################

def stacked_config(npsa, L, ll=1):
  '''
  Returns a stacked configuration of order npsa. The result is a list x of 
  size npsa, where x[i-1] is the position of the i-th particle. If npsa is 
  smaller or equal to Nmax, where Nmax is the maximum number of particles 
  that can fit onto the lattice, then the resulting configuration has npsa 
  particles, otherwise it has Nmax particles.
  '''
  
  Nmax = N_max(L, ll)
  if npsa > Nmax:
    npsa = Nmax
   
  # initialize list
  x = []
  
  if npsa > 0:
    # stack all particles
    i = 1
    xi = 1
    while (xi <= L and i <= npsa):
      x.append(xi)
      i += 1
      if i <= npsa:
        xi = 1 + (i-1)*ll
    
  return x
  


###############################################################################
#
# DEFINE: Returns the number of non-zero coefficients of a given PSA order.
#
###############################################################################

def nonzero_coeffs(n, L, ll=1):
  '''
  Returns the number of non-zero coefficients of a given order PSA order n, 
  given the lattice size L and the particle size ll.
  '''

  if n == 0:
    return 1

  elif n > 0:
    sum_i = np.zeros(n)  
    for ii in range(n):
      i = ii + 1
      term_i1 = choose(L-i*ll+i,i)
      sum_k = np.zeros(ll-1)    
      for kk in range(ll-1):
        k = kk + 1
        term_k = choose(L-k-(i-1)*ll+i-1,i-1)
        sum_k[kk] = term_k
      term_i2 = np.sum(sum_k)
      sum_i[ii] = term_i1 + term_i2

    Cn = 1 + np.sum(sum_i)
      
    return Cn


###############################################################################
#
# DEFINE: Returns the total number of all non-zero coefficients up to a given
# PSA order. 
#
###############################################################################

def total_coeffs(K, L, ll=1):
  '''
  Returns the total number of all non-zero coefficients up to a given PSA order
  K, given the lattice size L and the particle size ll. It is advised to use 
  this function if one is calling psa_compute with save_coeffs=True, in order to 
  check how much lines the resulting file will have.
  '''
  
  Nmax = N_max(L, ll)
  sum_Cn = np.zeros(K+1)
  for n in range(K+1):
    Cn = nonzero_coeffs(n, L, ll)
    sum_Cn[n] = Cn

  if K <= Nmax:

    result = np.sum(sum_Cn)
    return result.astype(int)

  elif K > Nmax:
  
    result = np.sum(sum_Cn) + ((K-Nmax)*nonzero_coeffs(Nmax, L, ll))
    return result.astype(int) 

  
  
###############################################################################
#
# DEFINE: Returns the next configuration in the PSA iteration.
#
###############################################################################

def next_config(xlist, L, ll=1):
  '''
  Returns next configuration in the PSA iteration from the input configuration 
  xlist. The input configuration is a list of size equal to the order of the 
  PSA. The elements of xlist are particle positions. The values of xlist for 
  particles that are missing are set to zero.
  
  For example, xlist = [1,2,0] means that the PSA order is 3, the 1st 
  particle is at lattice site 1, the 2nd particle is at lattice site 2, and 
  the 3rd particle is missing. On the other hand, xlist = [1,2] means that 
  the PSA order is 2, the 1st particle is at lattice site 1 and the 2nd 
  particle is at lattice site 2.
  '''
  
  # the order of the PSA, also the maximum number of particles for which 
  # the PSA coefficient is non-zero
  npsa = len(xlist)
  
  # number of zeros in xlist 
  nzeros = xlist.count(0)
  
  # number of particles in xlist
  N = npsa - nzeros
  
  if xlist[N-1] < L:
  
    # moves rightmost particle one lattice site to the right
    xlist[N-1] += 1
    
    # stackes npsa - N particles next to the rightmost particle, 
    # if possible
    if N < npsa:
      j = N
      xj = xlist[N-1] + ll
      while (xj <= L and j < npsa):
        xlist[j] = xj
        j += 1
        if j < npsa:
          xj = xlist[j-1]+ll
  else:
    xlist[N-1] = 0
  
  return xlist
  
  
  
###############################################################################
#
# DEFINE: Returns the total exit rate excluding initiation rate.
#
############################################################################### 
  
def e_0(xlist, wlist, ll=1): 
  '''
  Returns the total exit rate excluding initiation, see Eq. 30 in the main
  text.
  '''
  
  # number of particles
  N = len(xlist) - xlist.count(0)

  if N == 0:
    total = 0
  elif N == 1:
    pos_x = xlist[0]-1 # particle position shifted by -1
    total = wlist[pos_x]
  else:
    total = 0
    for particle_number in range(N-1):
      site = xlist[particle_number]
      next_site = xlist[particle_number+1] 
      if next_site-site > ll:
        total += wlist[site-1]
    site = xlist[N-1]
    total += wlist[site-1]
    
  return total
  
  
  
###############################################################################
#
# DEFINE: Computes the 1st order coefficient for a given configuration.
#
############################################################################### 
  
def c_1(xlist, wlist):
  '''
  Returns the 1st order coefficient c_1, see Eqs. 27 and 28 in the main text.
  '''

  if xlist[0] == 0: # empty configuration, implement Eq.(28)
    total = 0.0
    for w in wlist:
      total += 1.0/w
    coeff = - total

  else:  # non-empty configuration, implement Eq.(27)
    coeff=1.0/wlist[xlist[0]-1] #-1 since we start with position 0 
  
  return coeff
 
 
 
###############################################################################
#
# DEFINE: Performs the PSA up to a given order.
#
############################################################################### 
 
def psa_compute(wlist, K, ll=1, save_coeffs=False, coeffs_file='Pcoeff.csv'):
  '''
  Computes coefficients of the current and local density for orders 0,...,K.
  '''

  # finds lattice size
  L = len(wlist)

  # finds maximum number of particles
  Nmax = N_max(L, ll)
  
  # open file for storing probability coefficients 
  if save_coeffs == True:
    f=open(coeffs_file, 'w', newline='')
    writer = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
    writer.writerow([0,[],1.0]) # add zeroth order coefficient

  # initializes the dictionary
  c_n = {}

  # initializes the list for storing local density coefficients
  rhocoeff = []
  
  # initializes the list for storing current coefficients
  Jcoeff = []
  
  # order = 0
  Jcoeff.append(1.0) # the current coefficient is 1
  for site_i in range(L): 
    rhocoeff.append([0.0]) # the local density coefficients are 0

  # main loop over orders n = 1,...,K
  for n in range(1,K+1):
  
    c_n_minus1 = c_n.copy()
    c_n.clear()
    
    # resets current and local density coefficients
    Jcoeff_sum = 0
    rhocoeff_sum = [0 for site_i in range(L)]
    
    # stacked configuration
    X = stacked_config(n, L, ll)
    
    # number of particles in the stacked configuration
    N = len(X) 
    
    if n == 1: # order = 1
      
      # computes c_1 for the stacked configuration
      coeff = c_1(X, wlist)
    
    else: # order > 1
      xn = X[-1]-1 # position of the last particle shifted by -1
      if n <= Nmax:
        Xtemp = X[1:]  
      else:  
        Xtemp = X[1:]
        Xtemp.append(0)
      
      # computes c_n for the stacked configuration
      coeff = c_n_minus1[str(Xtemp)]/wlist[xn]
        
    # adds the stacked configuration coefficient to the dictionary c_n
    c_n[str(X)] = coeff
    
    # writes the stacked configuration coefficient to the file
    if save_coeffs == True:
      writer.writerow([n,np.trim_zeros(X,'b'),coeff])
          
    # updates the local density coefficients
    for particle_number in range(N):
      pos_x = particle_number*ll # particle position shifted by -1 
      rhocoeff_sum[pos_x] += coeff

    # iterates over all configurations belonging to the nth order
    while X[0] != 0: # iteration ends when the last particle leaves the lattice
    
      # selects next configuration
      X = next_config(X, L, ll)
      
      # number of particles in X
      N = len(X) - X.count(0)

      if n == 1: # order = 1
      
        # computes c_1 for configuration X
        coeff = c_1(X, wlist)
      
      else: # order > 1
      
        if X[0] == 0: # empty configuration
        
          # computes c_n for configuration X
          coeff = -sum(c_n.values())
          
        else:
          # total exit rate excluding initiation
          e_0X = e_0(X, wlist, ll)
          
          # computes term (a) in Eq.(29)
          a = 0
          if X[0] == 1 :
            if n <= Nmax:
              Xtemp = X[1:]  
            else:  
              Xtemp = X[1:]
              Xtemp.append(0)
          
            a = c_n_minus1[str(Xtemp)]

          # computes term (b1) in Eq.(29)
          b1 = 0
          if X[0] > 1 :
            Xtemp = X.copy()
            Xtemp[0] = Xtemp[0]-1
            pos_x1_min1 = Xtemp[0]-1 # -1 because wlist starts from index 0
            b1 = wlist[pos_x1_min1] * c_n[str(Xtemp)]
          
          # computes term (b2) in Eq.(29)
          b2 = 0
          if N > 1: # do this only if there is more than one particle
            for particle_number in range(N-1):
              Xtemp = X.copy()
              if X[particle_number + 1] - X[particle_number] > ll:
                Xtemp[particle_number+1] -= 1
                pos_x_min1 = Xtemp[particle_number+1]-1 #-1 again because wlist starts from index 0
                b2 +=  wlist[pos_x_min1] * c_n[str(Xtemp)]

          # computes term (c) in Eq.(29)
          c = 0
          if (X[N-1] <= L-ll) and (n >= N +1):
            Xtemp = X.copy()
            Xtemp[N] = L
            c = wlist[-1] * c_n[str(Xtemp)]

          # computes term (d) in Eq.(29)
          d = 0
          if (X[0] > ll) and (n-1 >= N):
            Xtemp = X.copy()
            if n <= Nmax:
              Xtemp.pop()
            d = c_n_minus1[str(Xtemp)]

          # computes c_n for configuration X
          coeff = (a + b1 + b2 + c - d)/e_0X
      
      # adds the computed coefficient to the dictionary c_n
      c_n[str(X)] = coeff
      
      # writes the coefficient to the file
      if save_coeffs == True:
        writer.writerow([n,np.trim_zeros(X,'b'),coeff])
      
      # updates the current coefficient
      if X[0] >= ll+1 or X[0] == 0:
        Jcoeff_sum += coeff
        
      # updates the local density coefficients
      for particle_number in range(N):
         pos_x = X[particle_number]-1 # particle position shifted by -1
         rhocoeff_sum[pos_x] += coeff
  
    Jcoeff.append(Jcoeff_sum)
 
    for site_i in range(L):
      rho_i_n = rhocoeff_sum[site_i]
      rhocoeff[site_i].append(rho_i_n)
      
  # closes the file
  if save_coeffs == True:
    f.close()
    
  return rhocoeff, Jcoeff
  
  
  
###############################################################################
#
# DEFINE: Returns local density for a given initiation rate. 
#
############################################################################### 

def local_density(rhocoeff, alpha):
  ''' 
  Returns local density up to every PSA order contained in rhocoeff. The
  output is a list of local density profiles rho such that rho[i] is the local
  density profile up to the PSA order i. The size of rho is the maximum PSA
  order plus 1.
  '''
  
  Kplus1 = len(rhocoeff[0]) # maximum order of the PSA increased by 1
  
  rho = [[] for order in range(Kplus1)]
  for coefficients in rhocoeff: # iteration over all lattice sites
    rho_sum = 0
    for order, coeff in zip(range(Kplus1), coefficients):
      rho_sum += alpha**order * coeff
      rho[order].append(rho_sum)
  
  return rho



###############################################################################
#
# DEFINE: Returns mean density given the local density.
#
############################################################################### 

def mean_density(local_density):

  mean_rho = []
  n = 0
  for local_density_n in local_density:
    rho_n = sum(local_density_n)/len(local_density_n)
    mean_rho.append(rho_n)
    n += 0
    
  return mean_rho



###############################################################################
#
# DEFINE: Returns particle current for a given initiation rate.
#
############################################################################### 

def current(Jcoeff, alpha):
  ''' 
  Returns particle current up to every PSA order contained in Jcoeff. The 
  output is a list of currents J such that J[i] is the values of the current
  up to the PSA order i. The size of J is the maximum PSA order plus 1.
  '''
  
  J = []
  J_sum = 0
  for order, coeff in zip(range(len(Jcoeff)), Jcoeff):
    J_sum += alpha**(order+1) * coeff
    J.append(J_sum)

  return J
