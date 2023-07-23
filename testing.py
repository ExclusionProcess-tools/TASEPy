import numpy as np
import matplotlib.pyplot as plt


from TASEPy import local_density
from TASEPy import mean_density
from TASEPy import current

from TASEPy import N_max
from TASEPy import stacked_config
from TASEPy import c_1
from TASEPy import next_config
from TASEPy import e_0


l = 10
K = 4
#wlist = [1.88, 1.52, 1.09, 1.38]*2
wlist = [7.6, 1.32, 5.65, 1.19, 1.59, 9.34, 6.74, 7.22, 6.42, 7.97, 1.42, 2.29, 3.18, 4.62, 2.7, 8.42, 2.89, 8.73, 7.3, 3.96, 9.92, 3.18, 2.61, 4.49, 9.74, 2.94, 9.64, 3.48, 4.51, 4.53, 4.85, 3.98, 2.31, 9.33, 1.0, 6.38, 8.15, 8.54, 5.77, 7.97, 4.09, 2.55, 9.72, 4.5, 3.76, 5.43, 7.46, 6.88, 2.67, 3.77, 1.42, 2.12, 2.03, 4.26, 4.39, 9.23, 9.27, 7.53, 1.24, 6.77, 7.21, 5.07, 3.78, 1.26, 9.97, 3.39, 4.24, 9.27, 8.73, 5.52, 1.44, 5.39, 5.08, 4.94, 5.14, 9.96, 5.07, 8.66, 2.66, 1.26, 2.34, 5.03, 4.24, 9.23, 3.49, 8.75, 5.43, 5.69, 3.57, 4.2, 2.38, 5.57, 3.82, 2.92, 8.14, 7.1, 8.66, 2.77, 8.24, 3.0, 5.98, 3.36, 3.6, 3.88, 4.87, 4.08, 7.25, 8.6, 6.8, 7.31, 9.57, 7.84, 5.94, 1.18, 9.36, 1.25, 7.22, 8.66, 7.3, 8.99, 1.4, 4.33, 2.77, 8.18, 6.7, 7.83, 2.91, 6.81, 5.43, 2.74, 4.37, 7.73, 6.21, 3.92, 6.17, 4.87, 7.96, 9.88, 2.75, 1.44, 1.25, 9.44, 2.65, 5.9, 6.25, 3.47, 1.14, 7.39, 8.09, 7.64, 7.01, 5.09, 1.16, 5.54, 2.21, 3.93, 9.0, 2.51, 6.87, 7.46, 8.88, 9.81, 4.65, 7.49, 1.24, 8.79, 2.68, 8.51, 9.19, 4.53, 1.87, 1.08, 8.64, 2.93, 6.18, 9.87, 7.92, 2.08, 1.08, 9.28, 4.19, 1.4, 5.18, 8.95, 5.97, 5.19, 2.37, 4.54, 4.94, 8.44, 8.92, 5.08, 7.85, 2.87, 7.39, 5.83, 8.94, 5.18, 1.02, 5.55, 1.04, 3.21, 2.3, 8.8, 6.82, 1.85, 4.69, 1.75, 5.58, 7.14, 6.51, 7.78, 2.96, 3.95, 6.79, 6.0, 5.71, 7.73, 6.89, 9.01, 3.44, 6.94, 8.39, 6.28, 2.88, 6.53, 8.95, 1.86, 2.11, 8.61, 4.39, 3.77, 6.17, 5.79, 4.17, 5.69, 8.07, 5.63, 3.53, 6.8, 1.1, 8.5, 2.54, 6.39, 1.46, 1.4, 1.83, 7.05, 6.33, 8.72, 7.19, 1.42, 4.19, 1.94, 1.4, 4.08, 9.51, 7.93, 7.03, 7.75, 3.94, 6.33, 9.05, 9.02, 9.28, 8.25, 7.57, 2.99, 6.89, 4.33, 1.81, 4.0, 8.18, 5.25, 4.36, 7.35, 6.28, 1.2, 3.29, 8.97, 8.39, 5.61, 7.41, 4.22, 9.43, 5.87, 3.62, 6.0, 7.99, 9.36, 6.37, 7.93, 2.27, 6.94, 4.71, 8.4, 8.19, 7.34, 6.8, 2.07]
L=len(wlist)

print(L)



def psa_compute(wlist, K, ll=1):
  '''
  Computes coefficients of the current and local density for orders 0,...,K.
  '''

  # finds lattice size
  L = len(wlist)

  # finds maximum number of particles
  Nmax = N_max(L, ll) 

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
          #
          


      
      # adds the computed coefficient to the dictionary c_n
      c_n[str(X)] = coeff
      
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
    
  return rhocoeff, Jcoeff




rhocoeff, Jcoeff = psa_compute(wlist, K, l)

alpha = 0.5
profile = local_density(rhocoeff, alpha)

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm


# Choose a colormap (you can change it to your preferred one)
color_map = cm.get_cmap('Blues')


# Create a new figure and axes
fig, ax = plt.subplots(figsize=(6, 4))

# Plot the data

# Define a list of linestyles
linestyles = ['-', '--', '-.', ':']

n = 0
for p in profile:

    line_color = color_map(0.2 + n / (K+1))
    # Cycle through the linestyles list
    linestyle = linestyles[n % len(linestyles)]

    if n != 0:
        ax.plot(profile[n], linewidth=2, label='n=' + str(n),
                color=line_color, linestyle=linestyle)
    n += 1

# Set the x and y axis labels
ax.set_xlabel(r'lattice site $i$', fontsize=12)
ax.set_ylabel(r'density $\rho_i$', fontsize=12)

# Set the title
ax.set_title(r'Density profile, $\alpha = $' +
             str(alpha), fontsize=14, fontweight='bold')

# Customize the tick parameters
ax.tick_params(axis='both', which='both', direction='in', labelsize=10)

# Customize the tick parameters for x-axis
# Set tick at every integer value
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
# Format tick labels as scalar values
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())


# Add a legend
ax.legend(fontsize=10)

# Add gridlines
ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

# Adjust the plot layout
plt.tight_layout()

# Save the plot as a high-resolution image
plt.savefig('figure_density.pdf', dpi=300)

# Show the plot
plt.show()
