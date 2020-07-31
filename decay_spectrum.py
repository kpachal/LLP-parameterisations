## Particle mass and energy must be set to determine other decay properties.
## Add angular distribution for more information.

## Should be callable with lists of particles too?


# Units in GeV
mass = 200

energy = 500

# Units in ns
taus = [0.01,0.1,1,10,100]

# Optional
# "isotropic", "central"
angular_dist = None

# And start
import numpy as np
from scipy import constants
import scipy.integrate as integrate
import matplotlib.pyplot as plt

for tau in taus :

# tau is in ns; convert to seconds
tau_s = tau * constants.nano
print("Tau in seconds:",tau_s)

## Mean lifetime tau:
# N(t) = N(0) exp(-t/tau)
## Mean distance xMean:
# N(x) = N(0) exp(-x/xMean)
# where xMean for relativistic particle = beta gamma c tau

# Energy is total energy; mass is rest mass in GeV/c^2
# Gamma = total energy/rest energy, where rest energy = rest mass * c^2
gamma = float(energy)/float(mass)
beta = np.sqrt(1.-(1./gamma)**2)
xMean = beta * gamma * constants.c * tau_s
print("Velocity is",beta,"of c.")
print("Beta-gamma is",beta*gamma)
print("ctau is",constants.c * tau_s,"m")
print("xMean is", xMean,"m")

## Now create decay distribution.
def N(x) :
  return np.exp(-x/xMean)

## Probability particle has decayed at location x
N_total = integrate.quad(N, 0, 200)[0]
def N_cumulative(x) :
  integral = integrate.quad(N,0,x)[0]
  return integral/N_total


## Eta distribution of particles
def eta_dist(theta) :
  # If we're not bothering with this effect, assuming everything is at eta=0
  if not angular_dist : return 0
#  elif angular_dist == "isotropic" :
    


#def eta_isotropic(theta) :
  

## Item 1: Plot number of particles versus decay distance, 
## and cumulative probability of decay.
# Reproduction of a Tova plot, as validation.
xvals = np.logspace(-6,2,200,endpoint=True)
yvals_cumulative = [N_cumulative(x) for x in xvals]
plt.figure()
plt.plot(xvals, N(xvals), 'b', xvals, yvals_cumulative, 'r')
plt.xlim(6e-6, 200)
plt.xscale('log')
plt.show()

